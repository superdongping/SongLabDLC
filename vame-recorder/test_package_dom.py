import os,subprocess,tempfile,time,re,json
from pathlib import Path
from urllib.request import urlopen,Request
exe=Path('dist/1.3.0/VAMERecorder/VAMERecorder.exe').resolve()
with tempfile.TemporaryDirectory(prefix='VAME_DOM_') as temp:
 p=subprocess.Popen([str(exe),'--port','43826','--data-dir',str(Path(temp)/'state'),'--synthetic','--no-browser'],creationflags=subprocess.CREATE_NO_WINDOW)
 token=''
 try:
  for _ in range(60):
   try:
    html=urlopen('http://127.0.0.1:43826',timeout=.5).read().decode();break
   except Exception:time.sleep(.2)
  token=re.search("const token='([^']+)'",html).group(1)
  env=os.environ.copy();env['VAME_TEST_OUTPUT']=temp
  r=subprocess.run(['node','test_dom.mjs'],env=env,capture_output=True,text=True,timeout=90)
  print(r.stdout,r.stderr);assert r.returncode==0
 finally:
  try:urlopen(Request('http://127.0.0.1:43826/api/shutdown',data=b'{}',headers={'X-VAME-Token':token}),timeout=10).read();p.wait(15)
  except Exception:p.terminate();p.wait()
print('PACKAGED_DOM_OK')
