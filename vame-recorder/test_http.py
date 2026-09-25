"""Black-box HTTP validation of the distributed executable, isolated test state."""
import argparse
import io
import json
from pathlib import Path
import re
import subprocess
import tempfile
import time
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import zipfile
import openpyxl

p=argparse.ArgumentParser();p.add_argument('--exe',default='dist/1.3.0/VAMERecorder/VAMERecorder.exe');a=p.parse_args()
root=Path(__file__).parent
temp=tempfile.TemporaryDirectory(prefix='VAME_HTTP_TEST_')
output=Path(temp.name)/'output';output.mkdir()
base='http://127.0.0.1:43825'
proc=subprocess.Popen([str(Path(a.exe).resolve()),'--port','43825','--data-dir',str(Path(temp.name)/'state'),'--synthetic','--no-browser'],creationflags=subprocess.CREATE_NO_WINDOW)
token=''
def req(route,data=None,auth=True):
    headers={'Content-Type':'application/json'}
    if auth:headers['X-VAME-Token']=token
    r=Request(base+route,headers=headers,data=None if data is None else json.dumps(data).encode())
    with urlopen(r,timeout=20) as f:return f.status,f.read()
def api(route,data=None):return json.loads(req('/api/'+route,data)[1])
def wait_status(sid,desired):
    deadline=time.monotonic()+25
    while time.monotonic()<deadline:
        s=next(s for s in api('state')['sessions'] if s['id']==sid)
        if s['status']==desired:return s
        if s['status'] in ('failed','interrupted'):raise AssertionError(s)
        time.sleep(.25)
    raise AssertionError('Timed out '+str(s))
try:
    for i in range(40):
        try:
            html=req('/',auth=False)[1].decode();break
        except URLError:time.sleep(.25)
    else:raise AssertionError('EXE failed to start')
    token=re.search("const token='([^']+)'",html).group(1)
    try:req('/api/state',auth=False);raise AssertionError('Missing auth accepted')
    except HTTPError as e:assert e.code==403
    assert api('state')['synthetic']
    api('output',{'path':str(output)})
    mouse=dict(mouse_id='HTTP_TEST',sex='M',cage='NO_ANIMAL',weight_g=25,group='Saline',operator='TEST')
    assert api('mouse',mouse)['result']['volume_ml']==.3125
    payload=dict(mouse_id='HTTP_TEST',phase='baseline',camera='TEST_SOURCE',size='640x480',fps=25,test=True,baseline_seconds=3,post_seconds=4)
    sid=api('start',payload)['result']
    # Reloading HTML does not stop the capture service.
    req('/',auth=False)
    s=wait_status(sid,'awaiting_injection')
    api('event',dict(session_id=sid,kind='injection',actual_volume_ml=.31,note='SIMULATED'))
    api('start',dict(session_id=sid,phase='post'))
    s=wait_status(sid,'completed')
    assert s['mouse']['group']=='Saline'
    assert s['actual_volume_ml']==.31
    assert s['phases']['post']['frames']==100
    assert all(q['mp4_status']=='queued' for q in s['phases'].values())
    api('queue',{'action':'process'})
    until=time.monotonic()+25
    while time.monotonic()<until:
        s=next(s for s in api('state')['sessions'] if s['id']==sid)
        if all(q['mp4_status']=='verified' for q in s['phases'].values()):break
        time.sleep(.1)
    assert all(q['mp4_status']=='verified' for q in s['phases'].values())
    raw=req('/api/export')[1]
    with zipfile.ZipFile(io.BytesIO(raw)) as z:assert z.testzip() is None
    workbook=openpyxl.load_workbook(io.BytesIO(raw),data_only=True)
    assert len(workbook.sheetnames)==13
    assert workbook['Cohort']['F8'].value==.3125
    assert workbook['Mouse 01']['B8'].value==.3125
    assert workbook['Mouse 01']['F8'].value in ('',None) # TEST is never exported as an actual experiment
    assert workbook['Mouse 01'].sheet_properties.pageSetUpPr.fitToPage
    api('shutdown',{})
    proc.wait(timeout=15)
    print('PACKAGED_HTTP_OK: auth, saving, timed phases, refresh, saline, XLSX cached formulas, shutdown')
finally:
    if proc.poll() is None:
        try:api('shutdown',{})
        except Exception:proc.terminate()
        proc.wait(timeout=15)
    temp.cleanup()
