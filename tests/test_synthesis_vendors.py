"""The Synthesis tab's company sections: IDT's client (idt.py) and where logins are kept (secret_store.py). No test reaches the network or
the real credential store: the HTTP call and the keyring are replaced."""
import base64
import io
import json
import sys
import urllib.error
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "app"))

import pytest

import idt
import secret_store
from user_error import UserError

CREDS = idt.Credentials("someone@lab.org", "p4ss-w0rd", "client-123", "s3cret-456")
SEQ = "ACGT" * 50   # 200 bases


def test_the_idt_login_and_screening_requests_have_the_shape_the_open_client_uses(monkeypatch):
    calls = []

    def post(url, data, headers, timeout=idt.TIMEOUT):
        calls.append((url, data, headers))
        return '{"access_token": "TOK", "expires_in": 3600}' if url == idt.TOKEN_URL else '[[{"Name": "openDNA", "Score": 4.5}, {"Name": "openDNA", "Score": 8}]]'

    monkeypatch.setattr(idt, "post", post)
    result = idt.screen(SEQ, CREDS)
    (url1, form, h1), (url2, body, h2) = calls
    assert url1 == idt.TOKEN_URL and h1["Authorization"] == "Basic " + base64.b64encode(b"client-123:s3cret-456").decode()
    assert dict(urllib.parse.parse_qsl(form.decode())) == {"grant_type": "password", "scope": "test", "username": "someone@lab.org", "password": "p4ss-w0rd"}
    assert url2 == idt.SCREEN_URL and h2["Authorization"] == "Bearer TOK" and json.loads(body) == [{"Name": "openDNA", "Sequence": SEQ}]
    assert result["score"] == 12.5 and len(result["rows"]) == 2 and json.loads(result["raw"])[0][1]["Score"] == 8


def test_an_idt_answer_is_read_whatever_its_nesting_and_a_missing_score_is_said():
    assert idt.parse('{"Score": 3}')["score"] == 3
    assert idt.parse('[[{"Score": 1}], [{"Score": 2.5}, {"Name": "x"}]]')["score"] == 3.5
    assert idt.parse('[{"Score": true}]')["score"] is None   # (a boolean is not a score)
    none = idt.parse('{"Message": "nothing to say"}')
    assert none["score"] is None and none["rows"] == [{"Message": "nothing to say"}]
    with pytest.raises(UserError, match="not JSON"):
        idt.parse("<html>maintenance</html>")


def test_a_length_idt_does_not_cover_is_refused_before_anything_is_sent(monkeypatch):
    monkeypatch.setattr(idt, "post", lambda *a, **k: pytest.fail("something was sent"))
    for n in (idt.MIN_BASES - 1, idt.MAX_BASES + 1):
        with pytest.raises(UserError, match=f"covers {idt.MIN_BASES} to {idt.MAX_BASES} bases.*Nothing was sent"):
            idt.screen("A" * n, CREDS)
    monkeypatch.setattr(idt, "post", lambda *a, **k: '{"access_token": "T"}' if a[0] == idt.TOKEN_URL else "[]")
    assert idt.screen("A" * idt.MIN_BASES, CREDS)["score"] is None and idt.screen("A" * idt.MAX_BASES, CREDS)["rows"] == []


def test_idt_errors_say_what_happened_without_repeating_the_login(monkeypatch):
    def refuse(request, timeout=None):
        raise urllib.error.HTTPError(request.full_url, 401, "Unauthorized", {}, io.BytesIO(b'{"error": "invalid_client"}'))

    monkeypatch.setattr(idt.urllib.request, "urlopen", refuse)
    with pytest.raises(UserError) as e:
        idt.screen(SEQ, CREDS)
    text = str(e.value)
    assert "401 Unauthorized" in text and "invalid_client" in text and "client secret" in text
    assert not any(secret in text for secret in (CREDS.password, CREDS.client_secret, CREDS.username, CREDS.client_id))

    def offline(request, timeout=None):
        raise urllib.error.URLError("no route")

    monkeypatch.setattr(idt.urllib.request, "urlopen", offline)
    with pytest.raises(UserError, match="Could not reach www.idtdna.com: no route"):
        idt.screen(SEQ, CREDS)
    monkeypatch.setattr(idt, "post", lambda *a, **k: '{"token_type": "Bearer"}')
    with pytest.raises(UserError, match="without an access token"):
        idt.screen(SEQ, CREDS)


class FakeKeyring:
    """A stand-in for the `keyring` module with a store that works (or fails to write)."""

    def __init__(self, module, fail_writes=False):
        backend = type("Backend", (), {"__module__": module})
        self.backend, self.store, self.fail_writes = backend(), {}, fail_writes

    def get_keyring(self):
        return self.backend

    def set_password(self, service, name, text):
        if self.fail_writes:
            raise secret_store.KeyringError("locked")
        self.store[(service, name)] = text

    def get_password(self, service, name):
        return self.store.get((service, name))

    def delete_password(self, service, name):
        self.store.pop((service, name), None)


def test_logins_go_to_the_systems_secure_store_or_stay_in_memory_never_in_a_file(monkeypatch):
    secret_store._session.clear()
    values = {"username": "u", "password": "pw"}
    good = FakeKeyring("keyring.backends.Windows")
    monkeypatch.setattr(secret_store, "keyring", good)
    assert secret_store.backend() == "Windows Credential Manager" and secret_store.secure()
    assert secret_store.save("idt", values) is True and json.loads(good.store[("openDNA", "idt")]) == values
    assert secret_store.load("idt") == values
    secret_store.forget("idt")
    assert secret_store.load("idt") == {} and good.store == {}

    for kr in (FakeKeyring("keyring.backends.fail"), None, FakeKeyring("keyring.backends.Windows", fail_writes=True)):
        monkeypatch.setattr(secret_store, "keyring", kr)
        assert secret_store.save("idt", values) is False   # (held until openDNA closes)
        assert secret_store.load("idt") == values
        secret_store.forget("idt")
        assert secret_store.load("idt") == {}
    monkeypatch.setattr(secret_store, "keyring", FakeKeyring("keyring.backends.fail"))
    assert secret_store.backend() is None and not secret_store.secure()
    secret_store._session.clear()
