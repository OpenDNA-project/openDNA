"""Reading UCSC .2bit genome files (what a track hub's `twoBitPath` names) and writing them as FASTA, the form the program opens.

No Qt. A .2bit file is a header, an index of (name, offset) and per sequence its length, its runs of N, its soft-masked (lower case)
runs and the bases packed four to a byte (T=0, C=1, A=2, G=3; the first base in the high bits). Version 1 (64-bit offsets, for
genomes over 4 GB) is read too. The FASTA is written in upper case: the soft-masking (repeats) is not kept.
"""
import os
import struct

from user_error import UserError

SIGNATURE = 0x1A412743
LINE = 80   # bases per FASTA line
_FOUR = ["".join("TCAG"[(b >> s) & 3] for s in (6, 4, 2, 0)).encode() for b in range(256)]   # a packed byte -> its four bases


def _header(f):
    """(byte order prefix, 64-bit offsets?, [(name, offset)]) of the open file."""
    head = f.read(16)
    if len(head) < 16:
        raise UserError("not a .2bit file (too short)")
    for order in ("<", ">"):
        if struct.unpack(order + "I", head[:4])[0] == SIGNATURE:
            break
    else:
        raise UserError("not a .2bit file (wrong signature)")
    version, count, _ = struct.unpack(order + "3I", head[4:])
    if version not in (0, 1):
        raise UserError(f".2bit version {version} is not supported (0 and 1 are)")
    index = []
    for _ in range(count):
        size = f.read(1)[0]
        name = f.read(size).decode("ascii", "replace")
        (offset,) = struct.unpack(order + ("Q" if version else "I"), f.read(8 if version else 4))
        index.append((name, offset))
    return order, index


def _blocks(f, order):
    (count,) = struct.unpack(order + "I", f.read(4))
    starts = struct.unpack(order + f"{count}I", f.read(4 * count))
    sizes = struct.unpack(order + f"{count}I", f.read(4 * count))
    return list(zip(starts, sizes))


def sequences(path):
    """[(name, length)] of the file, in its order."""
    out = []
    with open(path, "rb") as f:
        order, index = _header(f)
        for name, offset in index:
            f.seek(offset)
            out.append((name, struct.unpack(order + "I", f.read(4))[0]))
    return out


def read_sequence(f, order, offset):
    """The bases (upper case, N where the file has N) of the record at `offset`, as bytes."""
    f.seek(offset)
    (size,) = struct.unpack(order + "I", f.read(4))
    n_runs = _blocks(f, order)
    _blocks(f, order)   # the soft-mask runs: not kept
    f.read(4)           # reserved
    packed = f.read((size + 3) // 4)
    if len(packed) < (size + 3) // 4:
        raise UserError("the .2bit file is cut short")
    bases = bytearray(b"".join(map(_FOUR.__getitem__, packed))[:size])
    for start, length in n_runs:
        bases[start:start + length] = b"N" * length
    return bytes(bases)


def to_fasta(path, out, progress=None):
    """Write the .2bit file `path` as the FASTA file `out` (through a temporary file: a stopped or failed run leaves none).
    `progress(done, total, name)` is told after each sequence and may return False to stop (None is returned then).
    Returns [(name, length)]."""
    part = str(out) + ".part"
    done = []
    with open(path, "rb") as f, open(part, "wb") as w:
        order, index = _header(f)
        for i, (name, offset) in enumerate(index):
            bases = read_sequence(f, order, offset)
            w.write(b">" + name.encode("ascii", "replace") + b"\n")
            for at in range(0, len(bases), LINE):
                w.write(bases[at:at + LINE] + b"\n")
            done.append((name, len(bases)))
            if progress is not None and progress(i + 1, len(index), name) is False:
                w.close()
                os.remove(part)
                return None
    os.replace(part, out)
    return done
