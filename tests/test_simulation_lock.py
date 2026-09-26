"""The machine-wide lock that serialises MR-zero simulations (CPU only, no MR-zero)."""

import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest

from mess import _lock as sim

ROOT = Path(__file__).resolve().parents[1]
WAITING = "another MR-zero simulation is running"

# Holds the lock for a while in a separate process, then reports when it let go.
HOLDER = """
import sys, time
from mess._lock import simulation_lock
with simulation_lock(sys.argv[1], verbose=False):
    time.sleep(float(sys.argv[2]))
    released = time.time()
print(released, "MRzeroCore" in sys.modules)
"""


def test_acquire_and_release(tmp_path):
    assert sim.LOCK_FILE.parent == Path(tempfile.gettempdir())
    lock = tmp_path / "sim.lock"
    with sim.simulation_lock(lock):
        assert int(lock.read_text().split()[0]) == os.getpid()
    assert not lock.exists()
    with pytest.raises(RuntimeError):
        with sim.simulation_lock(lock):
            raise RuntimeError("simulation failed")
    assert not lock.exists()


def test_release_keeps_a_lock_taken_over(tmp_path):
    """A lock broken as stale and re-acquired elsewhere is not deleted by its old owner."""
    lock = tmp_path / "sim.lock"
    with sim.simulation_lock(lock):
        lock.write_text("12345 0123abcd")
    assert lock.read_text() == "12345 0123abcd"


def test_second_process_waits_for_release(tmp_path, capsys):
    lock = tmp_path / "sim.lock"
    holder = subprocess.Popen([sys.executable, "-c", HOLDER, str(lock), "1.5"], cwd=ROOT,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        deadline = time.time() + 120
        while not (lock.exists() and lock.read_text()):
            assert holder.poll() is None, holder.communicate()[1]
            assert time.time() < deadline, "the other process never took the lock"
            time.sleep(0.02)
        with sim.simulation_lock(lock, poll=0.05):
            acquired = time.time()
            assert int(lock.read_text().split()[0]) == os.getpid()
        out, err = holder.communicate(timeout=120)
    finally:
        holder.kill()
    assert holder.returncode == 0, err
    released, mrzero_imported = out.split()
    assert acquired >= float(released)
    assert mrzero_imported == "False"
    assert capsys.readouterr().out.count(WAITING) == 1
    assert not lock.exists()


@pytest.mark.parametrize("case", ["owner has exited", "too old", "empty and old"])
def test_stale_lock_is_taken_over(tmp_path, capsys, case):
    lock = tmp_path / "sim.lock"
    finished = subprocess.Popen([sys.executable, "-c", "pass"])
    finished.wait()             # keeping the Popen object keeps its PID from being reused
    content = {"owner has exited": f"{finished.pid} 0123abcd",
               "too old": f"{os.getpid()} 0123abcd", "empty and old": ""}[case]
    lock.write_text(content)
    if case != "owner has exited":
        old = time.time() - 3 * 3600
        os.utime(lock, (old, old))
    start = time.time()
    with sim.simulation_lock(lock, poll=0.05, stale_after=3600):
        assert int(lock.read_text().split()[0]) == os.getpid()
    assert time.time() - start < 5
    assert not lock.exists()
    assert not lock.with_name(lock.name + ".break").exists()
    assert WAITING not in capsys.readouterr().out


def test_pid_alive():
    finished = subprocess.Popen([sys.executable, "-c", "pass"])
    finished.wait()
    assert sim._pid_alive(os.getpid())
    assert not sim._pid_alive(finished.pid)
    assert not sim._pid_alive(0) and not sim._pid_alive(-1)
