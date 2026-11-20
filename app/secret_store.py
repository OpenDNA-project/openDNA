"""Passwords and API secrets, kept as safely as the system allows: the Windows Credential Manager, the macOS Keychain, or the Linux Secret
Service (GNOME Keyring, KWallet), through the `keyring` package. Never a file: where the system offers no secure store (a bare Linux, or
`keyring` missing) a secret lives in memory until openDNA closes, and `secure()` says so for the tab to tell the user."""
import json

try:
    import keyring
    from keyring.errors import KeyringError
except ImportError:   # (the package is a dependency, but a broken install must not break the window)
    keyring = None

    class KeyringError(Exception):
        pass

SERVICE = "openDNA"
_session = {}   # what could not be stored safely: this run only


def backend():
    """The name of the secure store in use ('Windows Credential Manager', ...), or None when there is none."""
    if keyring is None:
        return None
    try:
        kr = keyring.get_keyring()
    except Exception:
        return None
    module = type(kr).__module__
    if module.startswith("keyring.backends.fail") or ("chainer" in module and not getattr(kr, "backends", [])):
        return None
    return {"Windows": "Windows Credential Manager", "macOS": "macOS Keychain", "SecretService": "Secret Service", "kwallet": "KDE Wallet"}.get(
        module.rsplit(".", 1)[-1], type(kr).__name__)


def secure():
    return backend() is not None


def save(name, values):
    """Store the dict `values` under `name`. Returns True when it went to the secure store, False when it is only held until openDNA closes."""
    if backend() is not None:
        try:
            keyring.set_password(SERVICE, name, json.dumps(values))
            return True
        except (KeyringError, OSError):
            pass
    _session[name] = dict(values)
    return False


def load(name):
    """The dict stored under `name` (from the secure store, else from this run), or {}."""
    if backend() is not None:
        try:
            text = keyring.get_password(SERVICE, name)
        except (KeyringError, OSError):
            text = None
        if text:
            try:
                return json.loads(text)
            except ValueError:
                return {}
    return dict(_session.get(name, {}))


def forget(name):
    _session.pop(name, None)
    if backend() is not None:
        try:
            keyring.delete_password(SERVICE, name)
        except (KeyringError, OSError):
            pass
