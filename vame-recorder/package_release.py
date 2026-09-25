from pathlib import Path
import hashlib
import json
import shutil
import zipfile

root=Path(__file__).parent
dest=root/'dist/1.3.0/VAMERecorder'
for name in ('README.md','VALIDATION.md','CHANGELOG.md','THIRD_PARTY_NOTICES.txt','VAME_Experimental_Log.xlsx'):
    shutil.copy2(root/name,dest/name)
shutil.copytree(root/'licenses',dest/'licenses',dirs_exist_ok=True)
source=dest/'Source';source.mkdir(exist_ok=True)
for name in ('VAMERecorder.spec','build_windows.py','package_release.py','package.json','package-lock.json','app.py','background_jobs.py','test_idle.py','media_quality.py','test_media_quality.py','index.html','test_recorder.py','test_http.py','test_package_dom.py','test_relaunch.py','test_dom.mjs','validate_camera.py','VAME_Experimental_Log.xlsx'):
    shutil.copy2(root/name,source/name)
(source/'requirements.txt').write_text('imageio-ffmpeg==0.6.0\npyinstaller==6.22.3\nopenpyxl==3.1.5\n',encoding='utf-8')
files={str(p.relative_to(dest)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*') if p.is_file() and p.name!='SHA256.json'}
(dest/'SHA256.json').write_text(json.dumps(files,indent=2),encoding='utf-8')
target=root/'VAMERecorder_Windows_English_1.3.0.zip'
with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in dest.rglob('*'):
        if p.is_file():z.write(p,Path('VAMERecorder')/p.relative_to(dest))
with zipfile.ZipFile(target) as z:assert z.testzip() is None
print('PACKAGE_OK',target,round(target.stat().st_size/1e6,1),'MB')
