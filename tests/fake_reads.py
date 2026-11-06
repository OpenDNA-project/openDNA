"""Simulated sequencing of a plasmid with the usual cloning problems, to try the genome view.

    python tests/fake_reads.py [plasmid.dna|.gb] [out folder]
    (defaults: demo/pFG051_Chlor.gb, tests/fake_reads/)

Makes three "clones" of the plasmid, each a mix of variants (haplotypes) at given shares, and
sequences each one twice at ~30x:
- Element AVITI: 2x150 paired-end, ~350-bp inserts, very few errors;
- Oxford Nanopore: mostly full-length linearised plasmid reads, ~4 % errors, most of them
  one base too many or too few in homopolymers (the ONT signature).
Reads are mapped with the app's own recipes (tools.toml: bwa mem / minimap2 map-ont, in WSL on
Windows), so the BAMs are what Tools -> Map Reads would give. truth.tsv lists every change.
"""
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "app"))
import runner  # noqa: E402  (app module: recipes and the WSL bridge)
from opendna_core import Document  # noqa: E402

COMPLEMENT = str.maketrans("ACGTacgtN", "TGCAtgcaN")
DEPTH = 30


def rc(s):
    return s.translate(COMPLEMENT)[::-1]


# Events in reference coordinates, 1-based inclusive like GenBank (start, end):
#   snp (pos, base) | del (start, end) | ins (after pos, bases) | inv (start, end)
#   dup (start, end): tandem copy | is (after pos, length): foreign insert + 9-bp target duplication
def apply(ref, events):
    """The sequence of a haplotype: events applied right to left so positions stay valid."""
    seq = ref
    for kind, a, b, _ in sorted(events, key=lambda e: e[1], reverse=True):
        s = a - 1
        if kind == "snp":
            seq = seq[:s] + b + seq[s + 1:]
        elif kind == "del":
            seq = seq[:s] + seq[b:]
        elif kind == "ins":
            seq = seq[:a] + b + seq[a:]
        elif kind == "inv":
            seq = seq[:s] + rc(seq[s:b]) + seq[b:]
        elif kind == "dup":
            seq = seq[:b] + seq[s:b] + seq[b:]
        elif kind == "is":
            foreign = "".join(random.Random(b).choice("ACGT") for _ in range(b))
            seq = seq[:a + 9] + foreign + seq[a:]  # bases a+1..a+9 end up on both sides
    return seq


def other_base(ref, pos, rng):
    return rng.choice([c for c in "ACGT" if c != ref[pos - 1].upper()])


def homopolymer_deletion(ref, lo, hi):
    """1-bp deletion in the longest homopolymer between lo and hi (synthesis/frameshift error)."""
    best, run_start = (0, lo), lo
    for i in range(lo, hi):
        if ref[i].upper() != ref[i - 1].upper():
            run_start = i
        best = max(best, (i - run_start + 1, run_start))
    return ("del", best[1] + 1, best[1] + 1, f"1-bp deletion in a {best[0]}-base homopolymer (frameshift)")


def clones(ref, rng):
    """{clone name: [(share, [events])]} -- shares of each clone add up to 1."""
    snp_cds = ("snp", 3200, other_base(ref, 3200, rng), "point mutation in fosA (all reads)")
    snp_primer = ("snp", 960, other_base(ref, 960, rng), "mutation inside the oAC1346 primer (oligo synthesis error)")
    ins_promoter = ("ins", 1850, "T", "1-bp insertion in the pJExD promoter")
    del_rbs = ("del", 1900, 1908, "9-bp deletion at the gBlock/RBS junction")
    frameshift = homopolymer_deletion(ref, 1100, 1500)
    mixed_snp = ("snp", 1400, other_base(ref, 1400, rng), "mixed colony: mutation in half of the reads (eilR)")
    big_del = ("del", 2100, 2700, "601-bp deletion in 40 % of the plasmids (intein to terminator lost)")
    inversion = ("inv", 977, 1555, "eilR CDS inverted (fragment ligated backwards)")
    dup = ("dup", 2972, 3038, "tandem duplication of the J23115 promoter + RBS (Gibson misassembly)")
    insertion = ("is", 3300, 768, "768-bp foreign insertion in fosA in 25 % (IS-like transposon, 9-bp TSD)")
    return {
        "cloneA_point_errors": [(1.0, [snp_cds, snp_primer, ins_promoter, del_rbs, frameshift])],
        "cloneB_mixed_population": [
            (0.3, []), (0.3, [mixed_snp]), (0.2, [big_del]), (0.2, [big_del, mixed_snp]),
        ],
        "cloneC_rearrangements": [(0.75, [inversion, dup]), (0.25, [inversion, dup, insertion])],
    }


def pick(haplotypes, n):
    """Haplotype of read n: evenly spread (golden-ratio steps), so with only ~40 nanopore reads
    the shares still match truth.tsv instead of drifting by chance."""
    u, total = n * 0.6180339887 % 1, 0
    for share, h in haplotypes:
        total += share
        if u < total:
            return h
    return haplotypes[-1][1]


def window(seq, start, n):
    """n bases of a circular sequence from `start` (reads go through the origin)."""
    return (seq * (n // len(seq) + 2))[start:start + n]


def element(haplotypes, length, rng):
    """[(read1, read2, qual1, qual2)]: ~0.1 % substitutions, more towards the 3' end."""
    pairs = []
    for k in range(DEPTH * length // 300):
        h = pick(haplotypes, k)
        frag = window(h, rng.randrange(len(h)), max(160, min(len(h), int(rng.gauss(350, 50)))))
        if rng.random() < 0.5:
            frag = rc(frag)
        reads = []
        for r in (frag[:150], rc(frag)[:150]):
            bases, qual = list(r), []
            for i in range(len(bases)):
                if rng.random() < 0.0005 + 0.003 * i / 150:
                    bases[i] = rng.choice([c for c in "ACGT" if c != bases[i]])
                    qual.append("5")  # Q20
                else:
                    qual.append("I")  # Q40
            reads.append(("".join(bases), "".join(qual)))
        pairs.append((reads[0][0], reads[1][0], reads[0][1], reads[1][1]))
    return pairs


def nanopore(haplotypes, length, rng):
    """[(read, qual)]: 70 % full-length linearised plasmids, the rest broken; ~1 % substitutions,
    ~1 % random indels, and homopolymers of 4+ bases often one base off."""
    reads, total = [], 0
    while total < DEPTH * length:
        h = pick(haplotypes, len(reads))
        n = len(h) - rng.randrange(30) if rng.random() < 0.7 else rng.randrange(500, len(h))
        r = window(h, rng.randrange(len(h)), n)
        if rng.random() < 0.5:
            r = rc(r)
        out, qual, i = [], [], 0
        while i < len(r):
            j = i
            while j < len(r) and r[j] == r[i]:
                j += 1
            run = j - i
            if run >= 4 and rng.random() < 0.25:
                run += rng.choice((-1, 1))  # homopolymer miscalled
            for _ in range(run):
                x = rng.random()
                if x < 0.01:
                    out.append(rng.choice([c for c in "ACGT" if c != r[i]]))
                    qual.append(chr(33 + rng.randint(3, 9)))
                elif x < 0.015:
                    continue  # base lost
                else:
                    out.append(r[i])
                    qual.append(chr(33 + rng.randint(12, 30)))
                if rng.random() < 0.005:
                    out.append(rng.choice("ACGT"))  # extra base
                    qual.append(chr(33 + rng.randint(3, 9)))
            i = j
        reads.append(("".join(out), "".join(qual)))
        total += len(out)
    return reads


def map_reads(recipe, ref, reads1, reads2, out):
    files = {"ref": str(ref), "reads1": str(reads1), "reads2": str(reads2 or ""), "out": str(out)}
    steps = runner.steps_for(runner.recipes()[recipe], files, 4)
    done = subprocess.run(runner.command_line(runner.script(" && ".join(steps))), capture_output=True, text=True)
    if done.returncode != 0:
        sys.exit(f"{recipe} failed:\n{done.stdout}\n{done.stderr}")


def main():
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "demo" / "pFG051_Chlor.gb"
    folder = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "tests" / "fake_reads"
    folder.mkdir(parents=True, exist_ok=True)
    doc = Document(str(source))
    ref, name = doc.sequence().upper(), runner.safe_name(doc.name)
    fasta = folder / "ref.fa"
    fasta.write_text(f">{name}\n{ref}\n")
    rng = random.Random(55)
    truth = ["clone\tshare\tkind\tstart\tend\tchange\tdescription"]
    for clone, haplotypes in clones(ref, rng).items():
        for event in sorted({e for _, events in haplotypes for e in events}, key=lambda e: e[1]):
            share = sum(w for w, events in haplotypes if event in events)
            kind, a, b, text = event
            end = b if kind in ("del", "inv", "dup") else a
            change = {"snp": f"{ref[a - 1]}>{b}", "ins": f"+{b}", "is": f"+{b} bp"}.get(kind, f"{end - a + 1} bp")
            truth.append(f"{clone}\t{share:.0%}\t{kind}\t{a}\t{end}\t{change}\t{text}")
        seqs = [(w, apply(ref, events)) for w, events in haplotypes]
        r1, r2 = folder / f"{clone}_element_R1.fastq", folder / f"{clone}_element_R2.fastq"
        with open(r1, "w") as f1, open(r2, "w") as f2:
            for n, (a, b, qa, qb) in enumerate(element(seqs, len(ref), rng)):
                f1.write(f"@{clone}_e{n}/1\n{a}\n+\n{qa}\n")
                f2.write(f"@{clone}_e{n}/2\n{b}\n+\n{qb}\n")
        ont = folder / f"{clone}_ont.fastq"
        with open(ont, "w") as f:
            for n, (r, q) in enumerate(nanopore(seqs, len(ref), rng)):
                f.write(f"@{clone}_ont{n}\n{r}\n+\n{q}\n")
        map_reads("bwa_mem", fasta, r1, r2, folder / f"{clone}_element")
        map_reads("minimap2_ont", fasta, ont, None, folder / f"{clone}_ont")
        print(f"{clone}: {clone}_element.bam, {clone}_ont.bam")
    (folder / "truth.tsv").write_text("\n".join(truth) + "\n")
    print(f"done: {folder}")


if __name__ == "__main__":
    main()
