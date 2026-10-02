from pathlib import Path
import shutil,zipfile,hashlib,json
root=Path(__file__).parent;dest=root/'dist/1.2.5-r4/BehaviorHub'
for name in ['README.md','VALIDATION.md','THIRD_PARTY_NOTICES.txt']:shutil.copy2(root/name,dest/name)
shutil.copy2(root/'docs/Behavior_Hub_1.2.5-r4_Quick_Start.pdf',dest/'Behavior_Hub_1.2.5-r4_Quick_Start.pdf')
(dest/'licenses').mkdir(exist_ok=True)
for p in (root/'licenses').iterdir():
 if not p.name.startswith(('numpy','OpenCV')):shutil.copy2(p,dest/'licenses'/p.name)
src=dest/'Source';src.mkdir(exist_ok=True)
for name in ['app.py','project_store.py','background_jobs.py','test_projects.py','media_quality.py','index.html','presets.json','get_default_behavior_options.m','BehaviorHub.spec','test_behaviorhub.py','test_package.py','test_ui.mjs','package.json','package-lock.json','build_windows.py','package_release.py']:
 shutil.copy2(root/name,src/name)
shutil.copy2(root/'requirements.txt',src/'requirements.txt')
for name in ['video_output.py','test_video_output.py','focus_capture.py','focus_ui.js','test_focus.py','test_focus_hardware.py','test_focus_package.py','test_focus_ui.mjs']:
 shutil.copy2(root/name,src/name)
shutil.copytree(root/'native',src/'native',dirs_exist_ok=True)
shutil.copytree(root/'assets',src/'assets',dirs_exist_ok=True)
shutil.copytree(root/'docs',src/'docs',dirs_exist_ok=True)
files={str(p.relative_to(dest)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*') if p.is_file() and p.name!='SHA256.json'}
(dest/'SHA256.json').write_text(json.dumps(files,indent=2))
target=root/'BehaviorHub_Windows_English_1.2.5-r4.zip'
with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in dest.rglob('*'):
  if p.is_file():z.write(p,Path('BehaviorHub')/p.relative_to(dest))
with zipfile.ZipFile(target) as z:assert z.testzip() is None
print('PACKAGE_OK',target,round(target.stat().st_size/1e6,1),'MB')
