"""Explicit test-only real camera validation. No animal data or actual injection."""
import argparse
import json
from pathlib import Path
import time
from app import Recorder,run_ff

p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=15);a=p.parse_args()
root=Path(__file__).parent/'validation'
root.mkdir(exist_ok=True)
stamp=time.strftime('%Y%m%d_%H%M%S')
r=Recorder(root/('state_'+stamp))
r.set_output(str(root.resolve()))
r.mouse(dict(mouse_id='CAMERA_TEST',sex='Unknown',cage='NO_ANIMAL',group='KA',weight_g=25,operator='software validation',notes='No mouse; injection event is simulated.'))
sid=r.start(dict(phase='baseline',mouse_id='CAMERA_TEST',camera='HD Pro Webcam C920',size='1920x1080',fps=25,input_format='mjpeg',test=True,baseline_seconds=10,post_seconds=a.seconds))
def wait():
    while r.active:
        time.sleep(1)
wait()
s=r.session(sid)
assert s['status']=='awaiting_injection',s
r.event(dict(session_id=sid,kind='injection',actual_volume_ml=.3125,note='SIMULATED event for software test; no injection.'))
r.start(dict(phase='post',session_id=sid));wait()
assert s['status']=='completed',s
results={}
for phase in ('baseline','post'):
    f=Path(s['folder'])/(phase+'.mkv')
    result=run_ff(['-v','error','-i',str(f),'-map','0:v:0','-f','null','-'],timeout=max(90,a.seconds))
    results[phase]=dict(decode_returncode=result.returncode,decode_errors=result.stderr.decode('utf-8','replace'),**s['phases'][phase])
    assert result.returncode==0,results[phase]
report=dict(session=sid,seconds=a.seconds,results=results)
(root/('camera_'+stamp+'.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
