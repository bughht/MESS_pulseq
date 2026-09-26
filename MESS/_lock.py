"""A machine-wide lock that lets one MR-zero simulation run at a time.

:func:`mess.simulation.simulate` holds it while it simulates, so simulations
started from different processes or notebooks never run concurrently.
"""

import os
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

LOCK_FILE = Path(tempfile.gettempdir()) / "mess_mrzero.lock"


def _pid_alive(pid):
    """Whether process ``pid`` exists (on Windows ``os.kill(pid, 0)`` would kill it)."""
    if not 0 < pid < 2 ** 31:
        return False
    if os.name != "nt":
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:                                 # another user's process
            pass
        return True
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel32.GetExitCodeProcess.argtypes = (wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD))
    kernel32.CloseHandle.argtypes = (wintypes.HANDLE,)
    handle = kernel32.OpenProcess(0x1000, False, pid)           # PROCESS_QUERY_LIMITED_INFORMATION
    if not handle:
        return ctypes.get_last_error() == 5                     # ERROR_ACCESS_DENIED: it exists
    code = wintypes.DWORD()
    ok = kernel32.GetExitCodeProcess(handle, ctypes.byref(code))
    kernel32.CloseHandle(handle)
    return not ok or code.value == 259                          # STILL_ACTIVE


def _read_lock(path):
    """Content and age [s] of a lock file; ``(None, None)`` if there is none."""
    try:
        age = time.time() - os.stat(path).st_mtime
        with open(path) as f:
            return f.read(), age
    except FileNotFoundError:
        return None, None
    except PermissionError:                                     # Windows: being removed
        return "", 0.0


def _is_stale(text, age, stale_after, grace=10.0):
    """A lock is stale if its process has ended or it is older than ``stale_after`` [s]."""
    try:
        pid = int(text.split()[0])
    except (IndexError, ValueError):            # no PID yet (being written) or left by a crash
        return age > grace
    return age > stale_after or not _pid_alive(pid)


def _remove(path, content=None):
    """Delete a lock file (only if it still holds ``content``, when given).

    Retries for up to 5 s while another process has the file open (Windows).
    """
    for attempt in range(50):
        try:
            if content is not None:
                with open(path) as f:
                    if f.read() != content:
                        return
            os.remove(path)
            return
        except FileNotFoundError:
            return
        except PermissionError:
            if attempt == 49:
                raise
            time.sleep(0.1)


def _remove_stale(path, text, grace=10.0):
    """Delete the stale lock ``text`` unless it was replaced meanwhile.

    Only the holder of a second, short-lived lock file may do so, so that two
    waiters cannot both break the lock and one delete the other's new lock.
    Returns False if another process is breaking the lock.
    """
    breaker = path.with_name(path.name + ".break")
    try:
        os.close(os.open(breaker, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
    except OSError:
        _, age = _read_lock(breaker)
        if age is not None and age > grace:                     # left behind by a crash
            _remove(breaker)
        return False
    try:
        _remove(path, text)
    finally:
        _remove(breaker)
    return True


@contextmanager
def simulation_lock(path=None, poll=1.0, stale_after=2 * 3600, verbose=True):
    """Machine-wide lock that lets one MR-zero simulation run at a time.

    :func:`simulate` holds it while it simulates, so that simulations started
    from different processes or notebooks never share the GPU.  The lock is a file
    (``LOCK_FILE`` in the system temp directory unless ``path`` is given),
    created exclusively and holding the owner's PID; other processes wait,
    checking every ``poll`` seconds.  A lock whose process no longer exists,
    or that is older than ``stale_after`` seconds, is stale and taken over.
    Not re-entrant.
    """
    path = Path(LOCK_FILE if path is None else path)
    token = f"{os.getpid()} {os.urandom(4).hex()}"
    announced = False
    while True:
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            break
        except FileExistsError:
            pass
        except PermissionError:                                 # Windows: old lock being removed
            if os.name != "nt" or _read_lock(path)[0] is None:
                raise
        text, age = _read_lock(path)
        if text is None:                                        # released meanwhile
            continue
        stale = _is_stale(text, age, stale_after)
        if stale and _remove_stale(path, text):
            continue
        if verbose and not announced and not stale:
            owner = text.split()[0] if text.split() else "?"
            print(f"another MR-zero simulation is running (PID {owner}); waiting for {path}",
                  flush=True)
            announced = True
        time.sleep(poll)
    try:
        os.write(fd, token.encode())
    finally:
        os.close(fd)
    try:
        yield
    finally:
        _remove(path, token)                    # unless it was taken over as stale
