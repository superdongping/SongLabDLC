"""Opt-in packaged EXE/C920 API smoke test, isolated project and state."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import tempfile
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

def run(camera):
    exe=Path('dist/1.2.5-r4/BehaviorHub/BehaviorHub.exe').resolve()
    base='http://127.0.0.1:43837'
    token=''
    def call(route,data=None):
        req=Request(base+'/api/'+route,data=None if data is None else json.dumps(data).encode(),headers={'X-SongScope-Token':token})
        return json.load(urlopen(req,timeout=30))
    def launch(args):
        nonlocal token
        p=subprocess.Popen(args,creationflags=subprocess.CREATE_NO_WINDOW)
        end=time.monotonic()+15
        while time.monotonic()<end:
            try:
                html=urlopen(base,timeout=.5).read().decode()
                token=re.search("const token='([^']+)'",html).group(1)
                assert '1.2.5' in html and 'focusMode' in html
                return p
            except OSError:time.sleep(.1)
        p.terminate();p.wait();raise AssertionError('EXE failed to start')
    with tempfile.TemporaryDirectory(prefix='BehaviorHub_Focus_EXE_') as temp:
        args=[str(exe),'--port','43837','--data-dir',str(Path(temp)/'state'),'--no-browser']
        p=launch(args);original=None;c=None
        try:
            ds=[d for d in call('devices')['details'] if camera.lower() in d['name'].lower()]
            assert len(ds)==1,'Choose an unambiguous camera'
            d=ds[0];c=dict(camera=d['name'],camera_id=d['id'],size='1280x720',fps=25,input_format='mjpeg',focus={'mode':'auto'})
            original=call('focus',dict(c,action='read'))['result']
            call('project',dict(action='new',name='Packaged focus test',folder=temp))
            c['focus']={'mode':'manual','value':original['minimum']+original['step']}
            call('focus',dict(c,action='apply'));call('focus',dict(c,action='save'))
            sid=call('start',dict(c,assay='OFT',duration_seconds=2,test=True))['result']
            end=time.monotonic()+20
            while time.monotonic()<end:
                st=call('state')
                if not st['active']:break
                time.sleep(.1)
            s=next(s for s in st['sessions'] if s['id']==sid)
            assert s['status']=='completed',s.get('error')
            assert s['focus']['verified'] and s['focus']['flags']==2 and s['focus']['value']==c['focus']['value']
            project=st['project_file']
            call('shutdown',{});p.wait(10)
            p=launch(args);call('project',dict(action='open',path=project))
            st=call('state');assert st['focus_profiles'][c['camera_id']]['focus']==c['focus']
            assert not st['focus_status'],'Restart must not reuse stale verification'
            # Omitted focus request must still restore the saved manual profile.
            inherited=dict(c);inherited.pop('focus')
            call('preview',inherited);st=call('state')
            assert st['focus_status']['verified'] and st['focus_status']['value']==c['focus']['value']
            assert st['focus_status']['flags']==2
            call('stop-preview',{})
            bad=dict(c,focus={'mode':'manual','value':original['maximum']+original['step']})
            try:call('start',dict(bad,assay='OFT',duration_seconds=2,test=True))
            except HTTPError as exc:assert exc.code==400
            else:raise AssertionError('Focus failure did not block recording')
            st=call('state');assert not st['active']
            assert st['sessions'][-1]['status']=='failed'
            assert not list(Path(st['sessions'][-1]['folder']).glob('*.mkv'))
            print('PACKAGED_C920_FOCUS_OK: manual capture, save, restart, device profile restore, fresh validation, failure blocks video')
        finally:
            try:
                if original and c:
                    restore={'mode':'auto'} if original['flags']&1 else {'mode':'manual','value':original['value']}
                    result=call('focus',dict(c,focus=restore,action='apply'))['result']
                    assert result['verified']
                    call('stop-preview',{})
                    print('PACKAGED_TEST_ORIGINAL_MODE_RESTORED',restore['mode'])
            finally:
                try:call('shutdown',{});p.wait(10)
                except Exception:p.terminate();p.wait(5)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--camera',required=True)
    run(parser.parse_args().camera)
