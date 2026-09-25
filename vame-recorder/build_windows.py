from pathlib import Path
import os,shutil,subprocess,sys
import imageio_ffmpeg
root=Path(__file__).resolve().parent
if os.name!='nt':raise SystemExit('Build on Windows to produce the Windows executable.')
shutil.copy2(imageio_ffmpeg.get_ffmpeg_exe(),root/'ffmpeg.exe')
subprocess.run([sys.executable,'-m','PyInstaller','--noconfirm','--distpath','dist/1.3.0','--workpath','build/1.3.0','VAMERecorder.spec'],cwd=root,check=True)
subprocess.run([sys.executable,'package_release.py'],cwd=root,check=True)
