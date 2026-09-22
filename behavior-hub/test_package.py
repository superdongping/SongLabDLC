import subprocess,tempfile,time,json,re,os
from pathlib import Path
from urllib.request import urlopen,Request
exe=str(Path('dist/1.1.0/BehaviorHub/BehaviorHub.exe').resolve())
with tempfile.TemporaryDirectory(prefix='BehaviorHub_HTTP_') as temp:
 args=[exe,'--port','43836','--data-dir',str(Path(temp)/'state'),'--synthetic','--no-browser']
 p=subprocess.Popen(args,creationflags=subprocess.CREATE_NO_WINDOW)
 try:
  for i in range(60):
   try:
    with urlopen('http://127.0.0.1:43836',timeout=.5) as r:html=r.read().decode()
    break
   except Exception:time.sleep(.2)
  token=re.search("const token='([^']+)'",html).group(1)
  duplicate=subprocess.run(args,timeout=15);assert duplicate.returncode==0;assert p.poll() is None
  env=os.environ.copy();env['BEHAVIORHUB_TEST_OUTPUT']=temp
  result=subprocess.run(['node','test_ui.mjs'],env=env,capture_output=True,text=True,timeout=90)
  print(result.stdout,result.stderr);assert result.returncode==0
  def post(route,data):
   return json.load(urlopen(Request('http://127.0.0.1:43836/api/'+route,data=json.dumps(data).encode(),headers={'X-SongScope-Token':token})))
  post('calibration',dict(action='start',camera='synthetic',size='320x240',fps=25))
  time.sleep(1)
  frozen=post('calibration',{'action':'freeze'})['result'];assert frozen['width']>0
  post('calibration',dict(action='select',key='OFT',kind='rectangle',points=[[50,50],[800,50],[800,650],[50,650]]))
  post('start',dict(camera='synthetic',size='160x120',fps=25,assay='NPR',test=True,duration_seconds=1))
  with urlopen(Request('http://127.0.0.1:43836/api/alignment',headers={'X-SongScope-Token':token})) as r:assert not json.load(r)['enabled']
  time.sleep(2)
  print('PACKAGED_CALIBRATION_AND_RECORDING_PREEMPT_OK')
  for route in ['export','json']:
   with urlopen(Request('http://127.0.0.1:43836/api/'+route,headers={'X-SongScope-Token':token})) as r:
    assert len(r.read())>100
  with urlopen(Request('http://127.0.0.1:43836/api/state',headers={'X-SongScope-Token':token})) as r:project_file=json.load(r)['project_file']
  urlopen(Request('http://127.0.0.1:43836/api/shutdown',data=b'{}',headers={'X-SongScope-Token':token})).read();p.wait(10)
  p=subprocess.Popen(args,creationflags=subprocess.CREATE_NO_WINDOW)
  for i in range(60):
   try:
    with urlopen('http://127.0.0.1:43836',timeout=.5) as r:html=r.read().decode()
    break
   except Exception:time.sleep(.2)
  token=re.search("const token='([^']+)'",html).group(1)
  with urlopen(Request('http://127.0.0.1:43836/api/state',headers={'X-SongScope-Token':token})) as r:state=json.load(r)
  urlopen(Request('http://127.0.0.1:43836/api/project',data=json.dumps({'action':'open','path':project_file}).encode(),headers={'X-SongScope-Token':token})).read()
  with urlopen(Request('http://127.0.0.1:43836/api/state',headers={'X-SongScope-Token':token})) as r:state=json.load(r)
  assert state['sessions'][0]['status']=='completed';assert state['sessions'][0]['mouse']['mouse_id']=='OptionalMouse'
  print('PACKAGED_RESTART_AND_EXPORT_OK')
 finally:
  try:
   urlopen(Request('http://127.0.0.1:43836/api/shutdown',data=b'{}',headers={'X-SongScope-Token':token})).read();p.wait(10)
  except Exception:p.terminate();p.wait()
