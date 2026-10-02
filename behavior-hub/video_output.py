"""Project-relative MP4 destinations, fixed at recording start."""
from pathlib import Path
import re
from project_store import inside


def naming_mode(value):
    if value not in ('record_id', 'datetime_mouse', 'datetime_behavior'):
        raise ValueError('Choose a supported video naming option.')
    return value


def reserve_mp4(root, sessions, sid, timestamp, mouse_id, mode, test=False, behavior='', automatic_id='Auto_ID01'):
    naming_mode(mode)
    mouse = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', mouse_id).strip(' .')[:80]
    stem = sid if mode == 'record_id' else timestamp + ('_' + mouse if mouse else '')
    if mode == 'datetime_behavior':
        label = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', behavior).strip(' .')[:80]
        if not label:
            raise ValueError('Enter a behavioral test name.')
        stem = timestamp + '_' + label + '_' + (mouse or automatic_id)
    if test:
        stem = 'TEST_' + stem
    folder = inside(Path(root), 'MP4')
    folder.mkdir(exist_ok=True)
    used = {s.get('phases', {}).get('recording', {}).get('mp4_relative', '').casefold() for s in sessions}
    existing = {p.name.casefold() for p in folder.iterdir()}
    number = 1
    while True:
        name = stem + (f'_{number:03d}' if number > 1 else '') + '.mp4'
        relative = 'MP4/' + name
        if name.casefold() not in existing and relative.casefold() not in used:
            return relative
        number += 1


def mp4_destination(root, relative):
    path = inside(Path(root), relative)
    if path.parent != inside(Path(root), 'MP4') or path.suffix.lower() != '.mp4':
        raise ValueError('Invalid project MP4 destination.')
    return path
