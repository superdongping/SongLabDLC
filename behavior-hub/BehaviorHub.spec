# -*- mode: python ; coding: utf-8 -*-
import sys
from pathlib import Path
sys.path.insert(0,str(Path(SPECPATH)/"vendor"))


a = Analysis(
    ['app.py'],
    pathex=[str(__import__('pathlib').Path(SPECPATH)/'vendor')],
    binaries=[('ffmpeg.exe', '.')]+[(str(p),'numpy.libs') for p in (__import__('pathlib').Path(SPECPATH)/'vendor/numpy.libs').glob('*.dll')],
    datas=[('index.html', '.'), ('presets.json', '.')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='BehaviorHub',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='BehaviorHub',
)
