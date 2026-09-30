"""Opt-in short physical-camera smoke test; restores the starting focus mode.

Usage: python test_focus_hardware.py --camera C920
Uses disposable TEST recordings. Does not prove optical sharpness or USB recovery.
"""
import argparse
from pathlib import Path
import tempfile
import time
from app import Recorder, ROOT, device_details
from focus_capture import FocusBridge


def run(camera):
    matches=[d for d in device_details() if camera.lower() in d['name'].lower()]
    if len(matches)!=1: raise ValueError('Select one unambiguous physical camera name.')
    d=matches[0]
    c=dict(camera=d['name'],camera_id=d['id'],size='1280x720',fps=25,input_format='mjpeg')
    probe=FocusBridge(ROOT,c)
    original=dict(probe.report);probe.close()
    restore={'mode':'auto'} if original['flags'] & 1 else {'mode':'manual','value':original['value']}
    r=None
    try:
        with tempfile.TemporaryDirectory(prefix='BehaviorHub_Focus_Hardware_') as temp:
            r=Recorder(Path(temp)/'state')
            r.project_action(dict(action='new',folder=temp,name='Disposable C920 test'))
            value=min(original['maximum'],original['minimum']+original['step'])
            c['focus']={'mode':'manual','value':value}
            try:
                for i in range(3):
                    r.preview_start(c)
                    end=time.monotonic()+10
                    while not r.preview and time.monotonic()<end: time.sleep(.1)
                    assert r.preview, r.error
                    assert r.focus_status['verified'] and r.focus_status['value']==value
                    r.focus_action(dict(c,action='save'))
                    sid=r.start(dict(c,assay='OFT',duration_seconds=2,test=True))
                    end=time.monotonic()+20
                    while r.active and time.monotonic()<end: time.sleep(.1)
                    s=r.session(sid)
                    assert s['status']=='completed',s.get('error')
                    assert s['focus']['verified'] and s['focus']['flags']==2
                    assert s['focus']['value']==value
                    q=s['phases']['recording']
                    print('C920_MANUAL_PREVIEW_TO_RECORD_OK',i+1,'frames',q['frames'],'seconds',q['media_seconds'])
                bad=dict(c,focus={'mode':'manual','value':original['maximum']+original['step']})
                try: r.start(dict(bad,assay='OFT',duration_seconds=2,test=True))
                except ValueError as exc: print('C920_INVALID_FOCUS_BLOCKED',str(exc))
                else: raise AssertionError('Invalid focus unexpectedly allowed recording')
                assert not list(Path(r.db['sessions'][-1]['folder']).glob('*.mkv'))
                r.queue_action({'action':'process'})
                end=time.monotonic()+30
                while any(s['phases']['recording'].get('mp4_status') in ('queued','processing') for s in r.db['sessions']) and time.monotonic()<end: time.sleep(.1)
                for s in r.db['sessions'][:3]:
                    q=s['phases']['recording'];assert q['mp4_status']=='verified',q
                    print('C920_MP4_VERIFIED',q['decoded_frames'],q['observed_fps'],q['qc'])
            finally:
                r.request_stop('hardware test cleanup')
                if r.worker:r.worker.join(15)
                r.close_services()
    finally:
        restored=FocusBridge(ROOT,c,restore)
        try: print('ORIGINAL_FOCUS_RESTORED',restored.report['flags'],restored.report['value'])
        finally: restored.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--camera',required=True)
    run(p.parse_args().camera)
