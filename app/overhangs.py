"""Four-base Golden Gate overhangs: which sets are known to ligate faithfully, and a check for
the ones you pick.

Data: Potapov V, Ong JL, Kucera RB, et al. "Comprehensive Profiling of Four Base Overhang
Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly", ACS Synth Biol 2018, 7(11),
2665-2674, doi:10.1021/acssynbio.8b00333. The sets below are typed from the supplementary
Tables 3 and 4 of its bioRxiv version (doi:10.1101/322297): the junction overhangs of the
12-fragment (high and low fidelity) and 24-fragment lac cassettes that the authors built and
counted blue/white colonies for. The paper's full 256 x 256 ligation matrix is NOT included
(it is only in the supplementary spreadsheets); the check below is a rule of thumb built on
what the paper reports (mismatch ligation, most of all G:T pairs), not a measured fidelity.
"""

# 25 overhangs used for the high-fidelity 24-fragment cassette (Table S4); the 12-fragment
# high-fidelity set (Table S3) is a subset of it.
VALIDATED = ["GGAG", "GATA", "GGCA", "GGTC", "TCGC", "GAGG", "CAGT", "GTAA", "TCCA", "CACA", "GAAT", "ATAG",
             "AGTA", "ATCA", "TCTT", "AGGT", "CAAA", "AAGC", "GCAC", "CAAC", "AACG", "CGAA", "GTCT", "TCAG",
             "CCAT"]
# The paper's deliberately poor 12-fragment set (Table S3): kept as a test of the check below.
POOR_EXAMPLE = ["GGAG", "GGTC", "AGCA", "CAGT", "GGTA", "GAAT", "GGTT", "TCTT", "GGTG", "GCAC", "AGCG", "GTCT", "CCAT"]

_COMP = str.maketrans("ACGT", "TGCA")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}   # mismatches that can pair as G:T wobble


def rc(seq):
    return seq.translate(_COMP)[::-1]


def mismatches(a, b):
    """(count, wobble count) of the bases where overhang `a` differs from `b`. A mismatch is
    wobble-capable when it is a transition (A<->G, C<->T): the two then pair as G:T."""
    diff = [(x, y) for x, y in zip(a, b) if x != y]
    return len(diff), sum(p in TRANSITIONS for p in diff)


def pair_risk(a, b):
    """How likely overhangs a and b are to be joined to each other by mistake: 'bad' (equal or
    reverse complements: certain), 'warn' (one base apart) or None. Compared as written and as
    reverse complement (b could also pair with a's complement). Two overhangs from VALIDATED are
    exempt from the one-base rule: the paper built one assembly from all 25 of them."""
    risk = None
    for other in (b, rc(b)):
        m, _ = mismatches(a, other)
        if m == 0:
            return "bad"
        if m == 1 and not (a in VALIDATED and b in VALIDATED):
            risk = "warn"
    return risk


def assess(overhangs):
    """Checks for the overhangs of one assembly: [(name, level, text)] with level ok / warn / bad / info."""
    out, k = [], len(overhangs[0]) if overhangs else 4
    dup = sorted({o for o in overhangs if overhangs.count(o) > 1})
    out.append(("Overhangs unique", "bad" if dup else "ok", "repeated: " + ", ".join(dup) if dup else "all different"))
    pal = [o for o in overhangs if o == rc(o)]
    out.append(("Overhangs not palindromic", "bad" if pal else "ok",
                "a palindrome ligates to itself: " + ", ".join(pal) if pal else "none"))
    bad, warn = [], []
    for i, a in enumerate(overhangs):
        for b in overhangs[:i]:
            if a == b:
                continue   # reported as repeated above
            risk = pair_risk(a, b)
            if risk:
                (bad if risk == "bad" else warn).append(f"{b}/{a}")
    out.append(("Overhangs orthogonal", "bad" if bad else "warn" if warn else "ok",
                ("equal or reverse complements: " + ", ".join(bad) + "; " if bad else "")
                + ("one base apart, so a mismatched pair (often G:T) can ligate: " + ", ".join(warn) if warn else "")
                or "no two overhangs are equal, complementary or one base apart"))
    weak = [o for o in overhangs if len(set(o)) == 1 or set(o) <= set("GC") or set(o) <= set("AT")]
    out.append(("Overhang GC mix", "warn" if weak else "ok",
                "all G/C or all A/T (weak or sticky): " + ", ".join(weak) if weak else "mixed"))
    if k == 4:
        known = [o for o in overhangs if o in VALIDATED or rc(o) in VALIDATED]
        out.append(("Overhangs from the validated set", "ok" if len(known) == len(overhangs) else "info",
                    f"{len(known)} of {len(overhangs)} are among the 25 overhangs Potapov et al. 2018 built a "
                    "24-fragment assembly from" + ("" if len(known) == len(overhangs) else
                                                 "; not in it: " + ", ".join(o for o in overhangs if o not in known))))
    return out


def junction_checks(overhang, others):
    """The checks of one overhang among the `others` of the same assembly, for its junction's card:
    [(short name, level, text)]."""
    out = []
    equal = [q for q in others if q == overhang]
    out.append(("unique", "bad" if equal else "ok", "another junction has the same overhang" if equal else "no other junction has it"))
    out.append(("not palindromic", "bad" if overhang == rc(overhang) else "ok",
                "a palindrome ligates to itself" if overhang == rc(overhang) else "not its own reverse complement"))
    certain = [q for q in others if q != overhang and rc(q) == overhang]
    near = [q for q in others if q != overhang and rc(q) != overhang and pair_risk(overhang, q) == "warn"]
    out.append(("orthogonal", "bad" if certain else "warn" if near else "ok",
                f"reverse complement of {', '.join(certain)}: they will ligate" if certain else
                f"one base from {', '.join(near)}: a mismatched pair (often G:T) can ligate" if near else
                "no other overhang is equal, complementary or one base apart"))
    weak = len(set(overhang)) == 1 or set(overhang) <= set("GC") or set(overhang) <= set("AT")
    out.append(("GC mix", "warn" if weak else "ok", "all G/C or all A/T (weak or sticky)" if weak else "mixed"))
    if len(overhang) == 4:
        known = overhang in VALIDATED or rc(overhang) in VALIDATED
        out.append(("validated", "ok" if known else "info", "one of the 25 overhangs of Potapov et al. 2018" if known else
                    "not one of the 25 overhangs of Potapov et al. 2018 (no verdict, only a preference)"))
    return out


def junction_level(overhang, others):
    """ok / warn / bad for one junction's overhang among the `others` of the same assembly."""
    if overhang == rc(overhang) or any(pair_risk(overhang, q) == "bad" or overhang == q for q in others):
        return "bad"
    weak = len(set(overhang)) == 1 or set(overhang) <= set("GC") or set(overhang) <= set("AT")
    return "warn" if weak or any(pair_risk(overhang, q) for q in others) else "ok"


def preference(overhang, chosen):
    """Penalty for picking `overhang` after the ones in `chosen` (lower is better): used when
    openDNA chooses overhangs itself."""
    pen = 100 * (overhang == rc(overhang)) + 5 * (len(set(overhang)) == 1 or set(overhang) <= set("GC") or set(overhang) <= set("AT"))
    pen += 0 if (len(overhang) != 4 or overhang in VALIDATED) else 3
    for q in chosen:
        pen += {"bad": 100, "warn": 10}.get(pair_risk(overhang, q), 0) + 100 * (overhang == q)
    return pen


if __name__ == "__main__":
    assert not [c for c in assess(VALIDATED) if c[1] in ("bad", "warn")], assess(VALIDATED)
    assert [c for c in assess(POOR_EXAMPLE) if c[1] == "warn"], assess(POOR_EXAMPLE)
    print("ok")
