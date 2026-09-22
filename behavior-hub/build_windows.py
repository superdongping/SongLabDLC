"""Build the Windows package using the pinned dependencies; no user data is accessed."""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import imageio_ffmpeg

root = Path(__file__).resolve().parent
if os.name != 'nt':
    raise SystemExit('Build Behavior Hub on Windows to produce a Windows executable.')
shutil.copy2(imageio_ffmpeg.get_ffmpeg_exe(), root / 'ffmpeg.exe')
subprocess.run([sys.executable, '-m', 'PyInstaller', '--noconfirm', '--distpath', 'dist/1.1.1', '--workpath', 'build/1.1.1', 'BehaviorHub.spec'], cwd=root, check=True)
subprocess.run([sys.executable, 'package_release.py'], cwd=root, check=True)
