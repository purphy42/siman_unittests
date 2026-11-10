"""Viewer programs named in this project's project_conf.py."""

import ast
import os
import signal
import time
from pathlib import Path


SOURCE = Path(__file__).resolve().parent.parent / 'project_conf.py'
NAMES = ('PATH2JMOL', 'PATH2VESTA', 'PATH2OVITO')


def viewer_paths():
    """Read PATH2JMOL, PATH2VESTA, and PATH2OVITO without loading the project."""
    found = {}
    for line in SOURCE.read_text(encoding='utf-8').splitlines():
        code = line.split('#', 1)[0].strip()
        if '=' not in code:
            continue
        name, _, raw = code.partition('=')
        name = name.strip()
        if name not in NAMES:
            continue
        try:
            value = ast.literal_eval(raw.strip())
        except (SyntaxError, ValueError):
            continue
        if isinstance(value, str) and value:
            found[name] = value
    return found


def apply_viewer_paths():
    """Point SIMAN at the Li|LLZO viewer binaries."""
    from siman import header

    paths = viewer_paths()
    missing = [name for name in ('PATH2JMOL', 'PATH2VESTA') if name not in paths]
    if missing:
        raise FileNotFoundError(
            'project_conf.py does not set ' + ', '.join(missing))
    for name, value in paths.items():
        if not Path(value).is_file():
            raise FileNotFoundError(f'{name} does not exist: {value}')
        setattr(header, name, value)
    return paths


def viewer_pids():
    """Process ids of Jmol and VESTA, including processes they have started."""
    commands = {}
    parents = {}
    for entry in os.listdir('/proc'):
        if not entry.isdigit():
            continue
        pid = int(entry)
        try:
            command = Path(f'/proc/{pid}/cmdline').read_bytes().replace(b'\x00', b' ').decode('utf-8', 'replace')
            status = Path(f'/proc/{pid}/status').read_text(encoding='utf-8', errors='replace')
        except OSError:
            continue
        commands[pid] = command
        for line in status.splitlines():
            if line.startswith('PPid:'):
                parents[pid] = int(line.split()[1])
                break
    opened = {
        pid for pid, command in commands.items()
        if 'Jmol.jar' in command or '/VESTA' in command or 'VESTA-gtk' in command
    }
    changed = True
    while changed:
        changed = False
        for pid, parent in parents.items():
            if parent in opened and pid not in opened:
                opened.add(pid)
                changed = True
    return opened


def close_new_viewers(already_open, pause=1.2):
    """Close Jmol and VESTA processes that were not running before a test opened them."""
    time.sleep(pause)
    opened = viewer_pids() - set(already_open)
    for pid in opened:
        try:
            os.kill(pid, signal.SIGTERM)
        except OSError:
            pass
    time.sleep(0.4)
    for pid in opened:
        try:
            os.kill(pid, 0)
        except OSError:
            continue
        try:
            os.kill(pid, signal.SIGKILL)
        except OSError:
            pass
    return sorted(opened)
