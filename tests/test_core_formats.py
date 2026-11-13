"""File formats read by the Rust core: BAM reads / pileup / coverage / variants (against a hand count and a model on BAM
files written here with the standard library), GFF3 / GTF / BED tracks, GenBank feature locations, and indexed genome FASTA
files. Nothing here needs Qt.
"""
import random
import struct
import sys
import zlib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import opendna_core as core
from opendna_core import Document

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    return s.translate(COMP)[::-1]


def rand_dna(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


# ---- BAM ---------------------------------------------------------------------------------------------------------------


def _reg2bin(beg, end):
    """The BAM index bin of [beg, end) (SAM spec section 5.3)."""
    end -= 1
    for shift, offset in ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)):
        if beg >> shift == end >> shift:
            return offset + (beg >> shift)
    return 0


def _bgzf(data):
    """BGZF (blocked gzip) of `data`, with the empty end-of-file block the format asks for."""
    out = b""
    for i in range(0, len(data), 0xff00):
        chunk = data[i:i + 0xff00]
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = comp.compress(chunk) + comp.flush()
        size = len(body) + 25            # block size minus 1: 18 header + body + 8 trailer - 1
        out += struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, size)
        out += body + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    return out + bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


CIGAR_CODES = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
SEQ_CODES = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}


def read(name, pos, cigar, seq, flag=0, mapq=60, qual=None, ref_id=0, mate_ref=-1, mate_pos=-1, tlen=0):
    """One read: `cigar` like [(5, "S"), (45, "M")]; `qual` a list of base qualities, or None for 'not stored'."""
    return dict(name=name, pos=pos, cigar=cigar, seq=seq, flag=flag, mapq=mapq, qual=qual, ref_id=ref_id,
                mate_ref=mate_ref, mate_pos=mate_pos, tlen=tlen)


def write_bam(path, refs, reads):
    """A coordinate-sorted BAM written with nothing but the standard library: `refs` = [(name, length)]; reads sorted here."""
    text = b"@HD\tVN:1.6\tSO:coordinate\n"
    head = b"BAM\1" + struct.pack("<i", len(text)) + text + struct.pack("<i", len(refs))
    for name, length in refs:
        head += struct.pack("<i", len(name) + 1) + name.encode() + b"\0" + struct.pack("<i", length)
    body = b""
    for r in sorted(reads, key=lambda r: (r["ref_id"] if r["ref_id"] >= 0 else 1 << 30, r["pos"])):
        cigar = b"".join(struct.pack("<I", n << 4 | CIGAR_CODES[op]) for n, op in r["cigar"])
        ref_len = sum(n for n, op in r["cigar"] if op in "MDN=X")
        seq = r["seq"]
        packed = bytearray((len(seq) + 1) // 2)
        for i, base in enumerate(seq):
            packed[i // 2] |= SEQ_CODES[base] << (4 if i % 2 == 0 else 0)
        qual = bytes(r["qual"]) if r["qual"] is not None else b"\xff" * len(seq)
        name = r["name"].encode() + b"\0"
        end = r["pos"] + max(ref_len, 1)
        rec = struct.pack("<iiBBHHHiiii", r["ref_id"], r["pos"], len(name), r["mapq"], _reg2bin(r["pos"], end), len(r["cigar"]),
                          r["flag"], len(seq), r["mate_ref"], r["mate_pos"], r["tlen"])
        rec += name + cigar + bytes(packed) + qual
        body += struct.pack("<i", len(rec)) + rec
    Path(path).write_bytes(_bgzf(head + body))


def build_read(name, pos, ops, ref, subs=None, ins="", clip="T", **kw):
    """A read made from the reference along `ops` ("5S20M3D10M2I15M" as a list of (n, op)): substituted bases `subs`
    {reference position: base}, inserted bases `ins` (split over the I operations), soft-clipped bases `clip`."""
    subs, seq, r, ins_at = subs or {}, "", pos, 0
    for n, op in ops:
        if op == "M":
            seq += "".join(subs.get(r + k, ref[r + k]) for k in range(n))
            r += n
        elif op == "D":
            r += n
        elif op == "I":
            seq += ins[ins_at:ins_at + n]
            ins_at += n
        elif op == "S":
            seq += clip * n
    return read(name, pos, ops, seq, **kw)


def model_walk(r):
    """('M', ref position, base, quality) / ('I', ref position it goes before, bases) / ('D', ref position) for a read."""
    pos, q, events = r["pos"], 0, []
    for n, op in r["cigar"]:
        if op in "M=X":
            for k in range(n):
                qual = 255 if r["qual"] is None else r["qual"][q + k]
                events.append(("M", pos + k, r["seq"][q + k], qual))
            pos, q = pos + n, q + n
        elif op == "I":
            events.append(("I", pos, r["seq"][q:q + n]))
            q += n
        elif op in "DN":
            events += [("D", pos + k) for k in range(n)]
            pos += n
        elif op == "S":
            q += n
    return events


def usable(r):
    return not r["flag"] & (0x4 | 0x400 | 0x200)           # unmapped, duplicate, failed QC: ignored


def model_pileup(reads, ref_id, start, end, min_baseq=0, min_mapq=0):
    """Hand-counted pileup: per position [A, C, G, T, N, deleted, insertion before]."""
    counts = [[0] * 7 for _ in range(end - start)]
    for r in reads:
        if r["ref_id"] != ref_id or not usable(r) or r["mapq"] < min_mapq:
            continue
        for ev in model_walk(r):
            if not start <= ev[1] < end:
                continue
            slot = counts[ev[1] - start]
            if ev[0] == "M" and (ev[3] == 255 or ev[3] >= min_baseq):
                slot["ACGT".index(ev[2]) if ev[2] in "ACGT" else 4] += 1
            elif ev[0] == "D":
                slot[5] += 1
            elif ev[0] == "I":
                slot[6] += 1
    return counts


def model_depth(reads, ref_id, length):
    depth = [0] * length
    for r in reads:
        if r["ref_id"] == ref_id and usable(r):
            for ev in model_walk(r):
                if ev[0] == "M" and ev[1] < length:
                    depth[ev[1]] += 1
    return depth


@pytest.fixture(scope="module")
def ref_a():
    return rand_dna(random.Random(61), 2000)


@pytest.fixture(scope="module")
def hand_bam(tmp_path_factory, ref_a):
    """One read of every kind on chrA (2000 bp) and chrB (500 bp)."""
    ref = ref_a
    swap = lambda p, avoid=(): next(b for b in "ACGT" if b != ref[p] and b not in avoid)
    reads = [
        build_read("perfect", 10, [(50, "M")], ref, qual=[40] * 50),
        build_read("reverse", 30, [(50, "M")], ref, subs={45: swap(45)}, flag=0x10, qual=[40] * 50),
        build_read("softclip", 100, [(5, "S"), (45, "M")], ref, clip="T", qual=[30] * 50),
        build_read("deletion", 200, [(20, "M"), (5, "D"), (25, "M")], ref, qual=[30] * 45),
        build_read("insertion", 300, [(20, "M"), (3, "I"), (27, "M")], ref, ins="GGG", qual=[30] * 50),
        read("unmapped", 310, [], "ACGTACGTAC", flag=0x4, mapq=0, qual=[30] * 10),
        build_read("duplicate", 100, [(50, "M")], ref, flag=0x400, qual=[40] * 50),
        build_read("qcfail", 100, [(50, "M")], ref, flag=0x200, qual=[40] * 50),
        build_read("lowmapq", 500, [(50, "M")], ref, subs={510: swap(510)}, mapq=5, qual=[40] * 50),
        build_read("lowbq", 500, [(50, "M")], ref, subs={520: swap(520)}, qual=[40] * 20 + [10] + [40] * 29),
        build_read("hardclip", 700, [(10, "H"), (40, "M"), (10, "H")], ref, qual=[40] * 40),
        build_read("pair1", 900, [(50, "M")], ref, flag=0x1 | 0x40 | 0x20, mate_ref=1, mate_pos=100, tlen=0, qual=[40] * 50),
        build_read("pair2", 1000, [(50, "M")], ref, flag=0x1 | 0x40 | 0x20, mate_ref=0, mate_pos=1150, tlen=200, qual=[40] * 50),
        read("chrB_read", 10, [(30, "M")], "ACGT" * 7 + "AC", ref_id=1, qual=[40] * 30),
    ]
    path = tmp_path_factory.mktemp("bam") / "hand.bam"
    write_bam(path, [("chrA", 2000), ("chrB", 500)], reads)
    return core.Bam(str(path)), reads, path



@pytest.fixture(scope="module")
def random_bam(tmp_path_factory, ref_a):
    """300 reads of random position, length, strand, mapping quality, base qualities and CIGAR (with soft clips, deletions
    and insertions), with sequencing errors; plus a few that are ignored."""
    rng = random.Random(62)
    reads = []
    for i in range(300):
        pos = rng.randrange(0, 1900)
        ops = [(rng.randint(15, 60), "M")]
        if rng.random() < .15:
            ops += [(rng.randint(1, 6), "D"), (rng.randint(10, 30), "M")]
        if rng.random() < .15:
            ops += [(rng.randint(1, 4), "I"), (rng.randint(10, 30), "M")]
        if rng.random() < .3:
            ops = [(rng.randint(1, 8), "S")] + ops
        if rng.random() < .3:
            ops += [(rng.randint(1, 8), "S")]
        subs = {pos + rng.randrange(0, 60): rng.choice("ACGTN") for _ in range(rng.randint(0, 3))}
        ref_len = sum(n for n, op in ops if op in "MD")
        if pos + ref_len > 2000:
            continue
        n_bases = sum(n for n, op in ops if op in "MIS")
        flag = rng.choice([0, 0, 0x10, 0x10, 0x4, 0x400, 0x200]) if i % 10 == 0 else rng.choice([0, 0x10])
        if flag == 0x4:
            reads.append(read(f"r{i}", pos, [], rand_dna(rng, 20), flag=flag, mapq=0, qual=[30] * 20))
            continue
        r = build_read(f"r{i}", pos, ops, ref_a, subs={p: b for p, b in subs.items()}, ins=rand_dna(rng, 10), clip="C",
                       flag=flag, mapq=rng.choice([0, 10, 30, 60, 60, 60]),
                       qual=[rng.choice([2, 10, 20, 30, 40]) for _ in range(n_bases)])
        reads.append(r)
    path = tmp_path_factory.mktemp("bam") / "random.bam"
    write_bam(path, [("chrA", 2000)], reads)
    return core.Bam(str(path)), reads


def model_snps(reads, ref, min_depth, min_fraction, min_baseq=0, min_mapq=0):
    counts = model_pileup(reads, 0, 0, len(ref), min_baseq, min_mapq)
    out = []
    for pos, c in enumerate(counts):
        depth = sum(c[:6])
        if depth < max(min_depth, 1):
            continue
        want = "ACGT".index(ref[pos])
        n = max(c[s] for s in range(4) if s != want)
        if n > 0 and n >= min_fraction * depth:
            out.append((pos, ref[pos], n, depth))
    return out


@pytest.fixture(scope="module")
def snp_bam(tmp_path_factory, ref_a):
    """Ten reads over 100..175; at 120 five of them carry a SNP (50 %), at 125 one (10 %), at 130 three (30 %)."""
    reads = []
    alt = lambda p: next(b for b in "ACGT" if b != ref_a[p])
    for i in range(10):
        subs = {}
        if i < 5:
            subs[120] = alt(120)
        if i == 5:
            subs[125] = alt(125)
        if i in (6, 7, 8):
            subs[130] = alt(130)
        reads.append(build_read(f"s{i}", 100 + i // 2, [(70, "M")], ref_a, subs=subs, qual=[40] * 70, flag=0x10 if i % 2 else 0))
    path = tmp_path_factory.mktemp("bam") / "snp.bam"
    write_bam(path, [("chrA", 2000)], reads)
    return core.Bam(str(path)), reads, alt



def test_bam_reads_pileup_coverage_and_variants_match_a_model(hand_bam, random_bam, snp_bam, ref_a):
    """BAM with one read of every kind: sequences and the .bai; reads (soft clips, deletions, insertions, strand, ignored
    reads, pairing, windows); pileup and coverage against a hand count and a model on random reads; variant calling: a
    planted SNP and the thresholds, the reference base, deletions and insertions, random reads against the model."""
    # --- bam sequences and index written next to the file ---
    bam, _, path = hand_bam
    assert bam.sequences() == [("chrA", 2000), ("chrB", 500)]
    assert path.with_name(path.name + ".bai").exists(), "an unindexed sorted BAM gets its .bai written"
    with pytest.raises(OSError):
        bam.reads("chrZ", 0, 10, "", 10)
    # --- bam reads soft clips deletions insertions strand and ignored reads ---
    bam, _, _ = hand_bam
    got, more = bam.reads("chrA", 0, 2000, ref_a, 1000)
    assert not more
    by_name = {r[4]: r for r in got}
    # unmapped, duplicate and QC-failed reads are not listed; hard-clipped one is (it reads the aligned part)
    assert set(by_name) == {"perfect", "reverse", "softclip", "deletion", "insertion", "lowmapq", "lowbq", "hardclip", "pair1", "pair2"}
    p = by_name["perfect"]
    assert (p[0], p[1], p[2], p[3], p[5], p[6], p[7]) == (10, 60, 1, 60, [], [], [])
    r = by_name["reverse"]
    assert (r[0], r[1], r[2]) == (30, 80, -1) and r[8] & 0x10
    assert len(r[5]) == 1 and r[5][0][0] == 45 and r[5][0][1] != ref_a[45]
    s = by_name["softclip"]
    assert (s[0], s[1], s[5]) == (100, 145, []), "the 5 clipped T are not aligned: the read starts at 100 and has no mismatch"
    d = by_name["deletion"]
    assert (d[0], d[1], d[7]) == (200, 250, [(220, 225)]), "20 + 5 deleted + 25 = 50 reference bases"
    i = by_name["insertion"]
    assert (i[0], i[1], i[6], i[5]) == (300, 347, [(320, 3)], []), "inserted bases are not mismatches"
    h = by_name["hardclip"]
    assert (h[0], h[1]) == (700, 740)
    # pairing: flags, the mate position (0-based), the mate on another sequence, the insert size
    a, b = by_name["pair1"], by_name["pair2"]
    assert (a[9], a[10], a[11]) == (100, True, 0)
    assert (b[9], b[10], b[11]) == (1150, False, 200)
    assert by_name["perfect"][9:] == (None, False, 0)
    # a window gets only the reads that overlap it; limit sets the flag
    assert [r[4] for r in bam.reads("chrA", 100, 130, ref_a[100:130], 100)[0]] == ["softclip"]
    assert [r[4] for r in bam.reads("chrA", 0, 200, ref_a[:200], 2)[0]] == ["perfect", "reverse"]
    assert bam.reads("chrA", 0, 200, ref_a[:200], 2)[1] is True
    assert bam.reads("chrA", 0, 200, ref_a[:200], 4)[1] is False, "exactly `limit` reads: nothing more"
    assert [r[4] for r in bam.reads("chrB", 0, 100, "A" * 100, 10)[0]] == ["chrB_read"]
    # the reference is compared case-insensitively
    assert bam.reads("chrA", 0, 200, ref_a[:200].lower(), 100)[0][1][5] == bam.reads("chrA", 0, 200, ref_a[:200], 100)[0][1][5]
    # --- bam pileup hand count on the designed reads ---
    bam, reads, _ = hand_bam
    # 100..150: the soft-clipped read (45 aligned bases from 100) only; duplicate and QC-failed reads at the same place are ignored
    pile = bam.pileup("chrA", 100, 150)
    assert [sum(c[:4]) for c in pile] == [1] * 45 + [0] * 5
    assert all(c[4:] == [0, 0, 0] for c in pile)
    for i in range(45):
        assert pile[i]["ACGT".index(ref_a[100 + i])] == 1
    # deletion: 220..225 counted as deleted (and in the depth), insertion counted before 320
    d = bam.pileup("chrA", 215, 230)
    assert [c[5] for c in d] == [0] * 5 + [1] * 5 + [0] * 5
    assert [sum(c[:6]) for c in d] == [1] * 15
    ins = bam.pileup("chrA", 315, 325)
    assert [c[6] for c in ins] == [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    # everything, against the hand count
    assert bam.pileup("chrA", 0, 2000) == model_pileup(reads, 0, 0, 2000)
    # a window that runs past the end of the sequence is cut to it
    assert len(bam.pileup("chrA", 1990, 3000)) == 10
    # filters: the base with quality 10 is dropped at min_baseq 20, a mapq 5 read at min_mapq 20; the reads with no
    # filter-relevant values are untouched
    p0, pq, pm = bam.pileup("chrA", 500, 550), bam.pileup("chrA", 500, 550, min_baseq=20), bam.pileup("chrA", 500, 550, min_mapq=20)
    assert sum(sum(c[:4]) for c in p0) == 100
    assert sum(sum(c[:4]) for c in pq) == 99 and sum(p0[20][:4]) - sum(pq[20][:4]) == 1
    assert sum(sum(c[:4]) for c in pm) == 50 and pm == model_pileup(reads, 0, 500, 550, min_mapq=20)
    assert pq == model_pileup(reads, 0, 500, 550, min_baseq=20)
    assert sum(sum(c[:4]) for c in bam.pileup("chrA", 500, 550, min_baseq=41)) == 0
    assert sum(sum(c[:4]) for c in bam.pileup("chrA", 500, 550, min_mapq=61)) == 0
    # --- bam coverage bins equal a hand count ---
    bam, reads, _ = hand_bam
    depth = model_depth(reads, 0, 2000)
    # narrow bins (< 512 bases): counted exactly from the reads; mean depth of each bin
    for start, end, bins in ((0, 200, 20), (95, 155, 6), (200, 260, 60)):
        width = (end - start) // bins
        expect = [sum(depth[start + b * width:start + (b + 1) * width]) / width for b in range(bins)]
        assert bam.coverage("chrA", start, end, bins) == pytest.approx(expect, abs=1e-6), (start, end, bins)
    # wide bins (>= 4 blocks of 128): from the per-block totals computed once for the whole sequence
    assert not bam.has_depth("chrA")
    wide = bam.coverage("chrA", 0, 1024, 2)
    assert bam.has_depth("chrA")
    assert wide == pytest.approx([sum(depth[0:512]) / 512, sum(depth[512:1024]) / 512], abs=1e-6)
    assert bam.coverage("chrA", 0, 2048, 4) == pytest.approx([sum(depth[i * 512:(i + 1) * 512]) / 512 for i in range(4)], abs=1e-6)
    # an empty window or no bins
    assert bam.coverage("chrA", 50, 50, 3) == [0.0, 0.0, 0.0] and bam.coverage("chrA", 0, 100, 0) == []
    # --- bam pileup and coverage equal the model on random reads ---
    bam, reads = random_bam
    assert len(reads) > 250
    for min_baseq, min_mapq in ((0, 0), (20, 0), (0, 30), (30, 60), (41, 0)):
        assert bam.pileup("chrA", 0, 2000, min_baseq, min_mapq) == model_pileup(reads, 0, 0, 2000, min_baseq, min_mapq), (min_baseq, min_mapq)
        assert bam.pileup("chrA", 700, 913, min_baseq, min_mapq) == model_pileup(reads, 0, 700, 913, min_baseq, min_mapq)
    # keyword form too
    assert bam.pileup("chrA", 0, 100, min_baseq=20, min_mapq=10) == model_pileup(reads, 0, 0, 100, 20, 10)
    depth = model_depth(reads, 0, 2000)
    assert bam.coverage("chrA", 0, 2000, 100) == pytest.approx([sum(depth[i * 20:(i + 1) * 20]) / 20 for i in range(100)], abs=1e-6)
    assert bam.coverage("chrA", 0, 1024, 2) == pytest.approx([sum(depth[0:512]) / 512, sum(depth[512:1024]) / 512], abs=1e-6)
    # the depth the pileup reports (A..T, N, deleted) is not the coverage: deleted bases are in the first, not in the second
    assert any(c[5] for c in bam.pileup("chrA", 0, 2000))
    # --- variants planted snp at half is called and thresholds exclude it ---
    bam, reads, alt = snp_bam
    found = bam.variants("chrA", ref_a, min_depth=5, min_fraction=0.2)
    assert [v for v in found if v[0] in (120, 125, 130)] == [(120, "snp", ref_a[120], alt(120), 5, 10), (130, "snp", ref_a[130], alt(130), 3, 10)]
    assert not [v for v in found if v[0] == 125], "1 of 10 reads is 10 %: below the default 20 %"
    # the fraction: 50 % passes at 0.5, not at 0.6; 30 % passes at 0.3 and 0.2, not at 0.4
    at = lambda **kw: {v[0] for v in bam.variants("chrA", ref_a, **kw)}
    assert at(min_depth=5, min_fraction=0.5) == {120}
    assert at(min_depth=5, min_fraction=0.6) == set()
    assert at(min_depth=5, min_fraction=0.4) == {120}
    assert at(min_depth=5, min_fraction=0.05) == {120, 125, 130}
    # the depth: ten reads cover these bases, so 10 passes and 11 does not
    assert at(min_depth=10, min_fraction=0.5) == {120}
    assert at(min_depth=11, min_fraction=0.5) == set()
    # filters: with every base quality 40, min_baseq 41 leaves no depth at all
    assert bam.variants("chrA", ref_a, 1, 0.1, min_baseq=41) == []
    assert bam.variants("chrA", ref_a, 1, 0.1, min_mapq=61) == []
    # all of it equals the model, every position
    brief = [(p, r, n, d) for p, k, r, a, n, d in bam.variants("chrA", ref_a, min_depth=1, min_fraction=0.1) if k == "snp"]
    assert brief == model_snps(reads, ref_a, 1, 0.1)
    # --- variants reference base is compared and reported correctly ---
    bam, _, alt = snp_bam
    # Hand a reference that has the other base at 120: now the five reads that carry `alt` agree with it and the other five disagree
    swapped = ref_a[:120] + alt(120) + ref_a[121:]
    v = [x for x in bam.variants("chrA", swapped, min_depth=5, min_fraction=0.5) if x[0] == 120]
    assert v == [(120, "snp", alt(120), ref_a[120], 5, 10)]
    # the reference in lower case is the same reference
    assert bam.variants("chrA", ref_a.lower(), 5, 0.2) == bam.variants("chrA", ref_a, 5, 0.2)
    # a reference N at the position is no SNP at all, and a read N is not an alternative base
    with_n = ref_a[:120] + "N" + ref_a[121:]
    assert not [x for x in bam.variants("chrA", with_n, 5, 0.5) if x[0] == 120]
    with pytest.raises(OSError):
        bam.variants("chrZ", ref_a)
    # --- variants deletions and insertions and the random reads ---
    bam, _, _ = hand_bam
    found = bam.variants("chrA", ref_a, min_depth=1, min_fraction=0.5)
    got = {(v[0], v[1]): v for v in found}
    # the deletion read: 5 deleted bases merged into one variant (220..225), the reference bases listed
    assert got[(220, "del")] == (220, "del", ref_a[220:225], "", 1, 1)
    # the insertion read: GGG before 320
    assert got[(320, "ins")] == (320, "ins", "", "GGG", 1, 1)
    # the reverse read's substitution at 45: one of the two reads covering it (50 %)
    snp45 = [v for v in found if v[0] == 45]
    assert len(snp45) == 1 and snp45[0][1:3] == ("snp", ref_a[45]) and snp45[0][4:] == (1, 2)
    # the soft-clipped T bases are not variants: no call before 100 from that read
    assert not [v for v in found if 95 <= v[0] < 100]
    # Random reads against the model (SNPs only; the order of equal alternatives is not specified so the alt is checked by count)
    bam2, reads = random_bam
    snps = [v for v in bam2.variants("chrA", ref_a, 5, 0.3, 10, 20) if v[1] == "snp"]
    assert [(p, r, n, d) for p, k, r, a, n, d in snps] == model_snps(reads, ref_a, 5, 0.3, 10, 20)
    counts = model_pileup(reads, 0, 0, 2000, 10, 20)
    assert all(counts[p]["ACGT".index(a)] == n for p, k, r, a, n, d in snps)


@pytest.fixture(scope="module")
def track_of_many(tmp_path_factory):
    rng = random.Random(71)
    path = tmp_path_factory.mktemp("track") / "many.bed"
    rows = []
    for i in range(400):
        s = rng.randrange(0, 9_000)
        rows.append((s, s + rng.randint(1, 600), f"feat{i}"))
    path.write_text("".join(f"chr1\t{s}\t{e}\t{n}\n" for s, e, n in rows))
    return core.Track(str(path)), rows


def test_tracks_gff_bed_parsing_query_density_and_search(track_of_many, tmp_path):
    """GFF3 / GTF / BED tracks: escaped attributes, names, coordinates, 3 to 12 BED columns, headers, malformed lines refused;
    query equals a brute-force overlap (and the limit flag), density bins equal the overlap count, search is case-insensitive
    and limited."""
    # --- gff3 escaped attributes names and coordinates ---
    path = tmp_path / "t.gff3"
    path.write_text("##gff-version 3\n"
                    "chr1\tsrc\tregion\t1\t1000\t.\t+\t.\tID=chr1\n"                              # whole-sequence entries are skipped
                    "chr1\tsrc\tgene\t10\t100\t.\t+\t.\tID=g1;Name=big%2Cgene\n"                   # %2C = comma
                    "chr1\tsrc\tgene\t120\t180\t.\t-\t.\tID=g2;Name=a%3Bb%3Dc\n"                   # ; and = escaped
                    "chr1\tsrc\tgene\t200\t300\t.\t.\t.\tName=100%25 %C3%A9\n"                      # %25 = %, %C3%A9 = e acute
                    "chr1\tsrc\tCDS\t310\t320\t.\t+\t0\tID=c1;gene=araC\n"
                    "chr1\tsrc\texon\t330\t340\t.\t+\t.\tID=e1;gene_name=bad%ZZ;note=x\n"           # a bad escape stays as written
                    "chr1\tsrc\tmRNA\t350\t360\t.\t+\t.\tgene_id \"x1\"; transcript_id \"t1\"; Name \"gtfname\"\n"   # GTF quoting, Name first
                    "chr1\tsrc\ttRNA\t370\t380\t.\t+\t.\tcolor=red\n"                              # no name at all: the type
                    "chr1\tsrc\tgene\t400\t400\t.\t+\t.\tName=onebase\n"                            # one base long
                    "##FASTA\n>chr1\nACGT\n")
    t = core.Track(str(path))
    items, more = t.query("chr1", 0, 1000, 100)
    assert not more
    assert [i[3] for i in items] == ["big,gene", "a;b=c", "100% é", "araC", "bad%ZZ", "gtfname", "tRNA", "onebase"], "in start order"
    by_name = {i[3]: i for i in items}
    assert by_name["big,gene"][:3] == (9, 100, 1), "1-based inclusive 10..100 is 0-based 9..100"
    assert by_name["a;b=c"][:3] == (119, 180, -1)
    assert by_name["100% é"][2] == 0, "'.' strand is unknown"
    assert by_name["araC"][4] == "CDS" and by_name["araC"][:3] == (309, 320, 1), "no Name: the gene attribute names it"
    assert by_name["gtfname"][:3] == (349, 360, 1)
    assert by_name["tRNA"][4] == "tRNA"
    assert by_name["onebase"][:2] == (399, 400)
    assert "chr1" in [s for s, _ in t.sequences()] and t.sequences() == [("chr1", len(items))]
    assert not [i for i in items if i[4] in ("region", "chromosome", "source")]
    # Name beats gene_name, gene, locus_tag, gene_id, ID, transcript_id: the order of the preference list
    order = ["Name", "gene_name", "gene", "locus_tag", "gene_id", "ID", "transcript_id"]
    for k in range(len(order)):
        attrs = ";".join(f"{key}=v_{key}" for key in reversed(order[k:]))
        p = tmp_path / f"o{k}.gff3"
        p.write_text(f"c\ts\tgene\t1\t5\t.\t+\t.\t{attrs}\n")
        assert core.Track(str(p)).query("c", 0, 10)[0][0][3] == f"v_{order[k]}"
    # --- gff malformed lines are errors not crashes ---
    for i, line in enumerate(["chr1\ts\tgene\t10\t5\t.\t+\t.\tName=x",       # end before start
                              "chr1\ts\tgene\t0\t5\t.\t+\t.\tName=x",        # GFF is 1-based: 0 is no position
                              "chr1\ts\tgene\tten\t20\t.\t+\t.\tName=x",     # not a number
                              "chr1\ts\tgene\t1\t20\t.\t+",                  # too few columns
                              "just some text"]):
        p = tmp_path / f"bad{i}.gff3"
        p.write_text("##gff-version 3\n" + line + "\n")
        with pytest.raises(OSError):
            core.Track(str(p))
    with pytest.raises(OSError):
        core.Track(str(tmp_path / "no_such_file.gff3"))
    # --- bed with 3 to 12 columns ---
    full = ["chr1", "100", "200", "geneX", "500", "-", "110", "190", "255,0,0", "2", "10,20", "0,80"]
    for columns in range(3, 13):
        p = tmp_path / f"t{columns}.bed"
        p.write_text("\t".join(full[:columns]) + "\n" + " ".join(["chr2", "5", "15"] + full[3:columns]) + "\n")   # tabs, then spaces
        t = core.Track(str(p))
        (a,) = t.query("chr1", 0, 1000)[0]
        assert (a[0], a[1], a[4]) == (100, 200, "bed"), f"{columns} columns: BED is 0-based, end-exclusive: used as written"
        assert a[3] == ("geneX" if columns >= 4 else ""), f"{columns} columns: the name is column 4"
        assert a[2] == (-1 if columns >= 6 else 0), f"{columns} columns: the strand is column 6"
        (b,) = t.query("chr2", 0, 1000)[0]
        assert (b[0], b[1], b[2], b[3]) == (5, 15, a[2], a[3]), f"{columns} columns: space separated gives the same"
        assert t.sequences() == [("chr1", 1), ("chr2", 1)], columns
    # --- bed headers comments blank lines zero length and errors ---
    p = tmp_path / "t.bed"
    p.write_text("track name=x\nbrowser position chr1:1-100\n# comment\n\nchr1\t10\t10\tempty\nchr1\t20\t30\tn\t0\t.\nchr1\t40\t50\tp\t0\t+\n")
    t = core.Track(str(p))
    items = t.query("chr1", 0, 100)[0]
    assert [(i[0], i[1], i[2], i[3]) for i in items] == [(10, 10, 0, "empty"), (20, 30, 0, "n"), (40, 50, 1, "p")]
    # a zero-length feature is found in a window that contains its position, not in one that starts at it
    assert [i[3] for i in t.query("chr1", 5, 15)[0]] == ["empty"]
    assert t.query("chr1", 10, 20)[0] == []
    assert [i[3] for i in t.query("chr1", 0, 11)[0]] == ["empty"]
    for bad in ("chr1\tnot-a-number\t10\n", "chr1\t30\t10\n", "chr1\t30\n"):
        p.write_text(bad)
        with pytest.raises(OSError):
            core.Track(str(p))
    # --- track query equals a brute force overlap and the limit flag ---
    t, rows = track_of_many
    for start, end in ((0, 10_000), (1000, 1100), (4321, 4322), (8900, 9700), (20_000, 30_000)):
        expect = sorted(n for s, e, n in rows if s < end and e > start)
        items, more = t.query("chr1", start, end, 100_000)
        assert sorted(i[3] for i in items) == expect and not more
        assert [i[0] for i in items] == sorted(i[0] for i in items), "start order"
    everything = t.query("chr1", 0, 10_000, 100_000)[0]
    assert len(everything) == 400
    # limit: exactly `limit` items, and the flag says whether there were more
    items, more = t.query("chr1", 0, 10_000, 399)
    assert len(items) == 399 and more is True
    items, more = t.query("chr1", 0, 10_000, 400)
    assert len(items) == 400 and more is False
    assert t.query("chr1", 0, 10_000, 1)[1] is True
    assert t.query("chr2", 0, 10_000) == ([], False), "unknown sequence: nothing, no error"
    # default limit is 5000: the 400 come in one piece
    assert len(t.query("chr1", 0, 10_000)[0]) == 400
    # an empty Name= is no name: the gene attribute, then the ID, names the feature
    p = tmp_path / "empty_name.gff3"
    p.write_text("##gff-version 3\nc\ts\tgene\t10\t50\t.\t+\t.\tID=g1;Name=;gene=lacZ\n"
                 "c\ts\tgene\t60\t90\t.\t+\t.\tID=g2;Name=\n")
    assert [i[3] for i in core.Track(str(p)).query("c", 0, 100, 10)[0]] == ["lacZ", "g2"]
    # --- track density bins equal the overlap count and sum to the item count ---
    t, rows = track_of_many
    # Each bin counts the items overlapping it, so an item over several bins is counted in each of them
    for start, end, bins in ((0, 10_000, 100), (0, 10_000, 7), (2000, 6000, 40)):
        width = (end - start) / bins
        expect = [sum(1 for s, e, _ in rows if s < start + (b + 1) * width and e > start + b * width) for b in range(bins)]
        got = t.density("chr1", start, end, bins)
        assert got == expect, (start, end, bins)
        assert len(got) == bins
    # Items that each sit inside one bin: the bins add up to the number of items
    p = tmp_path / "tidy.bed"
    p.write_text("".join(f"c\t{s}\t{s + 5}\tf{s}\n" for s in range(2, 1000, 10)))      # 100 items of 5 bases in bins of 10
    tidy = core.Track(str(p))
    assert sum(tidy.density("c", 0, 1000, 100)) == 100
    assert tidy.density("c", 0, 1000, 100) == [1] * 100
    # No items in a window, no bins, an empty window
    assert tidy.density("c", 2000, 3000, 4) == [0, 0, 0, 0] and tidy.density("c", 0, 100, 0) == []
    assert tidy.density("c", 50, 50, 3) == [0, 0, 0]
    assert tidy.density("nope", 0, 100, 2) == [0, 0]
    # --- track search is case insensitive and limited ---
    p = tmp_path / "s.bed"
    p.write_text("chr1\t0\t10\tLacZ\nchr1\t20\t30\tlacY\nchr2\t0\t5\tÉCOLI-gene\nchr2\t10\t15\tother\n" + "".join(f"chr3\t{i}\t{i + 1}\tx{i}lac\n" for i in range(10)),
                 encoding="utf-8")
    t = core.Track(str(p))
    for needle in ("lac", "LAC", "LaC"):
        got = t.search(needle, 100)
        assert sorted(g[4] for g in got) == sorted(["LacZ", "lacY"] + [f"x{i}lac" for i in range(10)]), needle
    assert [g[4] for g in t.search("écoli", 10)] == ["ÉCOLI-gene"], "unicode case folding too"
    assert t.search("ecoli", 10) == [], "no accent folding"
    assert len(t.search("lac", 3)) == 3 and len(t.search("lac")) == 12
    row = t.search("lacz")[0]
    assert row == ("chr1", 0, 10, 0, "LacZ", "bed")
    assert t.search("zzz") == []


# ---- feature locations ---------------------------------------------------------------------------------------------------

# GenBank location text -> how it is written back (and the spans, 0-based end-exclusive, and the strand)
LOCATION_FORMS = [
    ("10..20", "10..20", [(9, 20)], 1),
    ("5", "5", [(4, 5)], 1),
    ("5..5", "5", [(4, 5)], 1),
    ("1..200", "1..200", [(0, 200)], 1),
    ("complement(10..20)", "complement(10..20)", [(9, 20)], -1),
    ("join(1..5,10..20)", "join(1..5,10..20)", [(0, 5), (9, 20)], 1),
    ("order(1..5,10..20)", "order(1..5,10..20)", [(0, 5), (9, 20)], 1),
    ("complement(join(1..5,10..20))", "complement(join(1..5,10..20))", [(0, 5), (9, 20)], -1),
    ("complement(order(1..5,10..20))", "complement(order(1..5,10..20))", [(0, 5), (9, 20)], -1),
    ("join(complement(1..5),complement(10..20))", "join(complement(1..5),complement(10..20))", [(0, 5), (9, 20)], -1),
    ("<10..20", "<10..20", [(9, 20)], 1),
    ("10..>20", "10..>20", [(9, 20)], 1),
    ("<1..>20", "<1..>20", [(0, 20)], 1),
    ("complement(<10..20)", "complement(<10..20)", [(9, 20)], -1),
    ("10^11", "10^11", [(10, 10)], 1),                                                     # the gap is where base 11 starts
    ("200^1", "200^1", [(0, 0)], 1),                                                       # across the origin of the 200 bp circle
    ("one-of(5,8)", "one-of(5,8)", [(4, 5), (7, 8)], 1),
    ("bond(1..5,10..20)", "bond(1..5,10..20)", [(0, 5), (9, 20)], 1),
    ("join(1..5,J00194.1:100..202)", "join(1..5,J00194.1:100..202)", [(0, 5)], 1),          # a part on another record
    ("join(1..5,gap(10),20..30)", "join(1..5,gap(10),20..30)", [(0, 5), (19, 30)], 1),      # a gap of 10 bases
]


def test_feature_locations_round_trip(tmp_path):
    """Feature locations: every GenBank form parses and is written back as a fixed point; bad ones are refused; in GenBank
    files a fuzzy end range is kept word for word and every form round trips through save."""
    # --- location parse and format round trip ---
    for text, written, spans, strand in LOCATION_FORMS:
        doc = Document.new_sequence("l", "ACGT" * 50, True)
        doc.add_features([("misc_feature", text, [("label", "x")])])
        f = doc.features()[-1]
        assert f.location == written and f.spans == spans and f.strand == strand, text
        # the written form is a fixed point: it parses to itself
        doc.add_features([("misc_feature", f.location, [("label", "y")])])
        g = doc.features()[-1]
        assert (g.location, g.spans, g.strand) == (f.location, f.spans, f.strand), text
        # set_feature takes the same text
        index = doc.set_feature(None, "gene", text, [])
        assert doc.features()[index].location == written, text
    # --- locations outside the sequence or not locations are refused ---
    doc = Document.new_sequence("l", "ACGT" * 50, True)
    for bad in ("10..300", "0..5", "20..10", "bad", ""):
        with pytest.raises((OSError, ValueError)):
            doc.add_features([("misc_feature", bad, [])])
    assert len(doc.features()) == 0, "a refused feature leaves nothing behind"
    with pytest.raises(OSError):
        doc.add_features([("two words", "1..5", [])])
    # --- a location with a fuzzy end range is kept word for word ---
    # `10..(20.30)` is not understood as a range, but it must never be cut to its first base: a file keeps it as written
    # (the feature has no extent, so it is not drawn), and typing one in is refused with a message.
    gb = ("LOCUS       t                        200 bp    DNA     circular UNK 01-JAN-1980\n"
          "FEATURES             Location/Qualifiers\n"
          "     misc_feature    10..(20.30)\n                     /label=\"f\"\n"
          "     misc_feature    join(1..5,(10.20)..(30.40))\n                     /label=\"g\"\n"
          "ORIGIN\n        1 " + "acgt" * 50 + "\n//\n")
    src = tmp_path / "fz.gb"
    src.write_text(gb)
    doc = Document(str(src))
    assert [f.location for f in doc.features()] == ["10..(20.30)", "join(1..5,(10.20)..(30.40))"]
    out = tmp_path / "fz2.gb"
    doc.write_genbank(str(out))
    assert [f.location for f in Document(str(out)).features()] == ["10..(20.30)", "join(1..5,(10.20)..(30.40))"]
    with pytest.raises(OSError, match="outside"):
        doc.add_features([("misc_feature", "10..(20.30)", [])])
    # --- locations from a genbank file round trip through save ---
    text = "".join(f"     misc_feature    {loc}\n                     /label=\"f{i}\"\n" for i, (loc, *_rest) in enumerate(LOCATION_FORMS)
                    if "gap" not in loc)
    gb = ("LOCUS       t                        200 bp    DNA     circular UNK 01-JAN-1980\nDEFINITION  .\nFEATURES             Location/Qualifiers\n"
          + text + "ORIGIN\n        1 " + "acgt" * 50 + "\n//\n")
    src = tmp_path / "t.gb"
    src.write_text(gb)
    doc = Document(str(src))
    before = [(f.location, f.spans, f.strand) for f in doc.features()]
    assert len(before) == len(LOCATION_FORMS) - 1
    out = tmp_path / "out.gb"
    doc.write_genbank(str(out))
    again = Document(str(out))
    assert [(f.location, f.spans, f.strand) for f in again.features()] == before
    # A location that spans two sides of the origin on a circle: join(190..200,1..10), with its strand
    doc.add_features([("misc_feature", "complement(join(190..200,1..10))", [])])
    f = doc.features()[-1]
    assert f.spans == [(189, 200), (0, 10)] and f.strand == -1


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    """Two chromosomes wrapped at 60 columns: one of 250 bases, one that is exactly two full lines (120)."""
    rng = random.Random(81)
    seqs = {"chr1": rand_dna(rng, 250), "chr2": rand_dna(rng, 120), "chr3": rand_dna(rng, 7)}
    path = tmp_path_factory.mktemp("genome") / "g.fa"
    path.write_text("".join(f">{n} some description\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n" for n, s in seqs.items()))
    return core.Fasta(str(path)), seqs, path


def test_genome_fasta_fetch_and_search_match_slicing_and_brute_force(genome, tmp_path):
    """Indexed genome FASTA: sequences and fetch (first / last base, line breaks, the end), random ranges equal slicing; the
    search across a line break, on both strands, with overlapping hits and against a brute-force scan; the IUPAC search on
    the file agrees with find_motif."""
    # --- genome fasta sequences and fetch first last base line breaks and the end ---
    fa, seqs, path = genome
    assert fa.sequences() == [("chr1", 250), ("chr2", 120), ("chr3", 7)]
    assert path.with_name(path.name + ".fai").exists()
    for name, s in seqs.items():
        n = len(s)
        assert fa.fetch(name, 0, 1) == s[0], "first base"
        assert fa.fetch(name, n - 1, n) == s[-1], "last base"
        assert fa.fetch(name, 0, n) == s, "everything"
        assert fa.fetch(name, 0, n + 100) == s, "an end past the sequence is cut to it"
        assert fa.fetch(name, n - 3, n + 100) == s[-3:]
    s = seqs["chr1"]
    # across line breaks: 59|60, 119|120, 179|180, a range of several lines
    for a, b in ((58, 62), (59, 61), (60, 61), (59, 60), (118, 122), (179, 181), (10, 245), (0, 60), (0, 61), (60, 120)):
        assert fa.fetch("chr1", a, b) == s[a:b], (a, b)
    # the last base of a sequence whose length is a multiple of the line length (no short last line)
    assert fa.fetch("chr2", 119, 120) == seqs["chr2"][119] and fa.fetch("chr2", 60, 120) == seqs["chr2"][60:]
    # empty or inverted ranges: nothing
    assert fa.fetch("chr1", 5, 5) == "" and fa.fetch("chr1", 9, 5) == "" and fa.fetch("chr1", 250, 260) == ""
    with pytest.raises(OSError):
        fa.fetch("chrZ", 0, 5)
    # --- genome fasta random ranges equal slicing ---
    fa, seqs, _ = genome
    rng = random.Random(82)
    for _ in range(300):
        name = rng.choice(list(seqs))
        a = rng.randrange(0, len(seqs[name]))
        b = rng.randrange(a, len(seqs[name]) + 5)
        assert fa.fetch(name, a, b) == seqs[name][a:b]
    # --- genome fasta search across a line break and motif search agrees with find motif ---
    fa, seqs, _ = genome
    s = seqs["chr1"]
    across = s[57:63]                                # six bases straddling the first line break
    assert any(p == 57 for p, st in fa.search("chr1", across) if st == 1)
    # the IUPAC search on the file equals the one on a string
    for motif, mism in (("GARTTC", 0), ("ACGT", 1), ("NNGG", 0)):
        on_file = fa.search_motif("chr1", motif, max_mismatches=mism)
        on_string = core.find_motif(s, motif, max_mismatches=mism)
        assert on_file == on_string
    # --- search on both strands, overlapping hits, against a brute-force scan ---
    rng = random.Random(83)
    filler = lambda n: rand_dna(rng, n, "AC")     # no G or T, so no accidental hit for a motif that has both
    motif = "GATTACA"
    seq = filler(30) + motif + filler(40) + rc(motif) + filler(25) + "AAAAAAAAA" + filler(5) + "GAATTC" + filler(10)
    path = tmp_path / "s.fa"
    path.write_text(">s\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    fa = core.Fasta(str(path))
    assert fa.search("s", motif) == [(30, 1), (77, -1)]
    assert fa.search("s", rc(motif)) == [(30, -1), (77, 1)], "the reverse complement of the motif finds the same places, strands swapped"
    assert fa.search("s", motif.lower()) == fa.search("s", motif), "case does not matter"
    # overlapping hits: AAAA in nine A is found at six positions, one after the other (and as TTTT on the other strand)
    a_start = seq.index("AAAAAAAAA")
    assert [p for p, st in fa.search("s", "AAAA") if a_start <= p <= a_start + 5] == list(range(a_start, a_start + 6))
    # a palindrome is reported once, on the top strand
    e = seq.index("GAATTC")
    assert fa.search("s", "GAATTC") == [(e, 1)]
    # against a brute-force scan, for several motifs, including ones across a line break (60)
    for m in ("AC", "ACA", "CA" * 5, "GAT", "TACAC"):
        expect = []
        for i in range(len(seq) - len(m) + 1):
            w = seq[i:i + len(m)]
            if w == m:
                expect.append((i, 1))
            elif w == rc(m):
                expect.append((i, -1))
        assert fa.search("s", m) == expect, m
    assert fa.search("s", "") == [] and fa.search("s", "G" * 500) == []
