from pathlib import Path
import shutil,zipfile,hashlib,json
root=Path(__file__).parent;dest=root/'dist/BehaviorHub'
for name in ['README.md','VALIDATION.md','THIRD_PARTY_NOTICES.txt']:shutil.copy2(root/name,dest/name)
shutil.copytree(root/'licenses',dest/'licenses',dirs_exist_ok=True)
src=dest/'Source';src.mkdir(exist_ok=True)
for name in ['app.py','media_quality.py','index.html','presets.json','get_default_behavior_options.m','BehaviorHub.spec','test_behaviorhub.py','test_package.py','test_ui.mjs','package.json','package-lock.json','build_windows.py','package_release.py']:
 shutil.copy2(root/name,src/name)
(src/'requirements.txt').write_text('imageio-ffmpeg==0.6.0\npyinstaller==6.22.3\n')
files={str(p.relative_to(dest)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*') if p.is_file() and p.name!='SHA256.json'}
(dest/'SHA256.json').write_text(json.dumps(files,indent=2))
target=root/'BehaviorHub_Windows_English_1.0.1.zip'
with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in dest.rglob('*'):
  if p.is_file():z.write(p,Path('BehaviorHub')/p.relative_to(dest))
with zipfile.ZipFile(target) as z:assert z.testzip() is None
print('PACKAGE_OK',target,round(target.stat().st_size/1e6,1),'MB')
