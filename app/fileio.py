"""Writing a file the safe way: write a temporary file, then `replace` it over the real one (a crash never leaves half a file).

Every such rename in the app goes through `replace`. On Windows a rename over a file is refused (PermissionError) while another
thread still has that file open, and the app has such readers: the background search for primer sites reads the other records of
a file, a genome search reads the reference. The reader lets go in a moment, so the rename is retried for a few seconds
(core/src/document.rs does the same for the saves it makes)."""
import os
import time
from pathlib import Path


def replace(tmp, path, tries=50, pause=0.1):
    """os.replace(tmp, path), retried for up to ~5 s when the target is a file that is busy; any other error is raised at once."""
    for _ in range(tries - 1):
        try:
            return os.replace(tmp, path)
        except PermissionError:
            target = Path(path)
            if not target.is_file() or not os.access(target, os.W_OK):   # a folder, or read-only: waiting will not help
                break
            time.sleep(pause)
    return os.replace(tmp, path)
