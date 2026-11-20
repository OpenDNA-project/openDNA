"""IDT's complexity screen for gBlocks gene fragments: would IDT make this sequence? Needs the user's own IDT API account (a username and
password, and the client id and secret IDT gives for API use). The sequence and these credentials go to idtdna.com, only when the user
presses the button (synthesis_panel.py asks first).

!!! NOT VERIFIED AGAINST THE LIVE SERVICE. IDT's API documentation is behind an IDT login and could not be read when this was written.
The addresses, the form of the request and the shape of the answer below follow what the open-source SBOL-utilities client
(`calculate_complexity_scores`, https://sbol-utilities.readthedocs.io/) is reported to do: a password-grant token from the identity server
(client id and secret as HTTP Basic), then the screening call with the sequences as a list of {"Name", "Sequence"}, answered with a list of
complexity findings per sequence, each with a "Score" (the sequence's score is the sum of its findings'). Check TOKEN_URL, SCREEN_URL and
`parse` against IDT's own documentation and one real call (the tab shows the raw reply to make that easy) before relying on it. !!!"""
import base64
import json
import urllib.error
import urllib.parse
import urllib.request
from typing import NamedTuple

import structure
from user_error import UserError

TOKEN_URL = "https://www.idtdna.com/Identityserver/connect/token"
SCREEN_URL = "https://www.idtdna.com/restapi/v1/Complexities/ScreenGblockSequences"
SIGN_UP_URL = "https://www.idtdna.com/pages/tools/apidoc"
MIN_BASES, MAX_BASES = 125, 3000   # the length range the gBlock screen is reported to cover (a score of 0 is returned outside it)
# (product, highest score IDT is quoted as accepting): gBlocks under 10, genes in custom vectors under 19. Sources differ: read IDT's page.
PRODUCTS = (("gBlocks Gene Fragment", 10.0), ("Gene in a custom vector", 19.0))
TIMEOUT = 120


class Credentials(NamedTuple):
    username: str
    password: str
    client_id: str
    client_secret: str


def check_length(seq):
    """UserError, before anything is sent, when IDT's screen does not cover this many bases."""
    if not MIN_BASES <= len(seq) <= MAX_BASES:
        raise UserError(f"IDT's screen covers {MIN_BASES} to {MAX_BASES} bases (this is {len(seq):,}): it would answer 0, which means nothing. Nothing was sent.")


def post(url, data, headers, timeout=TIMEOUT):
    """The text IDT answers to a POST. UserError says what went wrong; it never repeats what was sent."""
    request = urllib.request.Request(url, data=data, headers={"User-Agent": structure.USER_AGENT, **headers})
    host = url.split("/")[2]
    try:
        with urllib.request.urlopen(request, timeout=timeout) as reply:
            return reply.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace").strip()[:300]
        hint = " Check the username, password, client id and client secret." if e.code in (400, 401, 403) else ""
        raise UserError(f"{host} answered {e.code} {e.reason}." + (f" Its message: {body}" if body else "") + hint) from None
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        raise UserError(f"Could not reach {host}: {getattr(e, 'reason', e)}") from None


def token(creds):
    basic = base64.b64encode(f"{creds.client_id}:{creds.client_secret}".encode()).decode()
    form = urllib.parse.urlencode({"grant_type": "password", "scope": "test", "username": creds.username, "password": creds.password}).encode()
    reply = post(TOKEN_URL, form, {"Authorization": f"Basic {basic}", "Content-Type": "application/x-www-form-urlencoded"})
    try:
        return json.loads(reply)["access_token"]
    except (ValueError, KeyError, TypeError):
        raise UserError("IDT answered the login without an access token (see IDT's API documentation: the login call may have changed).") from None


def findings(answer):
    """The dictionaries found in the answer, whatever its nesting (a list per sequence of lists of findings)."""
    if isinstance(answer, dict):
        return [answer]
    if isinstance(answer, list):
        return [row for part in answer for row in findings(part)]
    return []


def parse(text):
    """{'raw': the answer as indented text, 'rows': its findings, 'score': sum of their "Score" numbers or None when there are none}."""
    try:
        answer = json.loads(text)
    except ValueError:
        raise UserError("IDT's answer is not JSON: " + text.strip()[:200]) from None
    rows = findings(answer)
    scores = [r["Score"] for r in rows if isinstance(r.get("Score"), (int, float)) and not isinstance(r.get("Score"), bool)]
    return {"raw": json.dumps(answer, indent=2), "rows": rows, "score": sum(scores) if scores else None}


def screen(seq, creds):
    """Ask IDT about `seq` (bases only): log in, send it, parse the answer. UserError for anything the user can fix or should know."""
    check_length(seq)
    bearer = token(creds)
    body = json.dumps([{"Name": "openDNA", "Sequence": seq}]).encode()
    return parse(post(SCREEN_URL, body, {"Authorization": f"Bearer {bearer}", "Content-Type": "application/json"}))
