"""Shared checks for the text a user types or pastes: each returns why a value is refused, in a sentence, or "" when it is fine,
so every field explains itself the same way (see docs/INPUT_AUDIT.md)."""
import re

API_KEY = re.compile(r"[A-Za-z0-9]*")


def text_problem(text, what, max_len, one_line=True):
    """A name or label: printable characters only (a pasted tab, line break or control character would end up in a file),
    at most `max_len`."""
    if one_line and any(not c.isprintable() for c in text):
        bad = sorted({repr(c)[1:-1] for c in text if not c.isprintable()})
        return f"The {what} has a tab, line break or control character ({' '.join(bad[:4])}): use plain letters, digits and spaces."
    if len(text) > max_len:
        return f"The {what} is too long: at most {max_len} characters ({len(text):,} given)."
    return ""


def api_key_problem(key):
    """An NCBI API key: letters and digits only (NCBI gives 36 hexadecimal characters)."""
    if not API_KEY.fullmatch(key):
        bad = sorted({c for c in key if not c.isalnum() or not c.isascii()})
        shown = " ".join(repr(c)[1:-1] if not c.isprintable() or c == " " else c for c in bad[:6])
        return f"An NCBI API key has only letters and digits; this one also has: {shown}. Copy it again from your NCBI account."
    if len(key) > 100:
        return "An NCBI API key is 36 characters; this one is far longer."
    return ""
