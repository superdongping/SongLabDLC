"""Regression tests for the Windows executable launcher; uses no real camera."""
import hashlib
import json
from pathlib import Path
import re
import socket
import subprocess
import tempfile
import time
from urllib.request import urlopen, Request

EXE=Path('dist/1.3.0/VAMERecorder/VAMERecorder.exe').resolve()
processes=[]
def launch(folder,port):
    p=subprocess.Popen([str(EXE),'--data-dir',str(folder),'--port',str(port),'--no-browser','--synthetic'],creationflags=subprocess.CREATE_NO_WINDOW)
    processes.append(p)
    return p
def ready(port):
    end=time.monotonic()+15
    while time.monotonic()<end:
        try:
            with urlopen(f'http://127.0.0.1:{port}/health',timeout=.5) as r:
                if json.load(r).get('version')=='1.3.0':return
        except Exception:pass
        time.sleep(.1)
    raise AssertionError('Service did not become ready')
def api(port,route,data=None):
    base=f'http://127.0.0.1:{port}'
    html=urlopen(base,timeout=3).read().decode()
    token=re.search("const token='([^']+)'",html).group(1)
    req=Request(base+'/api/'+route,headers={'X-VAME-Token':token,'Content-Type':'application/json'},data=json.dumps(data).encode() if data is not None else None)
    with urlopen(req,timeout=5) as r:return json.load(r)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
try:
    with tempfile.TemporaryDirectory(prefix='VAME_RELAUNCH_') as t:
        folder=Path(t)/'state';output=Path(t)/'videos';output.mkdir()
        first=launch(folder,43827);ready(43827)
        api(43827,'mouse',dict(mouse_id='RELAUNCH_TEST',cage='TEST',sex='F',group='KA',weight_g=25,notes='User data: 原始备注'))
        before=digest(folder/'registry.json')
        # Simulate closing the tab: leave the service alive, then launch again.
        for _ in range(3):
            duplicate=launch(folder,43827)
            assert duplicate.wait(10)==0
            assert first.poll() is None
            assert digest(folder/'registry.json')==before
        # Same data directory, different requested port must reopen the real port.
        duplicate=launch(folder,43828);assert duplicate.wait(12)==0
        assert digest(folder/'registry.json')==before
        # Relaunch during capture must not recover or interrupt the live session.
        api(43827,'output',{'path':str(output)})
        sid=api(43827,'start',dict(mouse_id='RELAUNCH_TEST',phase='baseline',test=True,camera='TEST',size='640x480',fps=25,baseline_seconds=6,post_seconds=2))['result']
        duplicate=launch(folder,43827);assert duplicate.wait(10)==0
        end=time.monotonic()+12
        while time.monotonic()<end:
            s=api(43827,'state')['sessions'][-1]
            if s['status']=='awaiting_injection':break
            assert s['status'] not in ('interrupted','failed'),s
            time.sleep(.2)
        assert s['status']=='awaiting_injection'
        assert s['phases']['baseline']['frames']==150
        api(43827,'event',dict(session_id=sid,kind='abandon',note='Test finished'))
        saved=digest(folder/'registry.json')
        api(43827,'shutdown',{});assert first.wait(12)==0
        # Full service shutdown followed by immediate restart.
        restarted=launch(folder,43827);ready(43827)
        assert digest(folder/'registry.json')==saved
        assert api(43827,'state')['mice'][0]['notes']=='User data: 原始备注'
        api(43827,'shutdown',{});assert restarted.wait(12)==0
        # Cold-start race: one server, one clean handoff.
        a=launch(folder,43827);b=launch(folder,43827);ready(43827)
        end=time.monotonic()+12
        while a.poll() is None and b.poll() is None and time.monotonic()<end:time.sleep(.1)
        assert sorted([a.poll() is None,b.poll() is None])==[False,True]
        assert next(p for p in (a,b) if p.poll() is not None).returncode==0
        api(43827,'shutdown',{})
        assert a.wait(12)==0 and b.wait(12)==0
        # An unrelated listener must not be treated as VAME or overwritten.
        with socket.socket() as sock:
            sock.bind(('127.0.0.1',43829));sock.listen()
            conflict=launch(Path(t)/'conflict',43829)
            assert conflict.wait(12)==1
        assert not (Path(t)/'conflict/registry.json').exists()
    print('RELAUNCH_OK: repeated launch, different port, during recording, shutdown/restart, cold-start race, unrelated port, data preservation')
finally:
    for p in processes:
        if p.poll() is None:
            p.terminate();p.wait(10)
