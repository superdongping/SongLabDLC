import io
import json
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import zipfile
import xml.etree.ElementTree as E
from app import Recorder, volume, export_workbook

class RecorderTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)
        self.r=Recorder(self.root/'state',synthetic=True)
        self.r.set_output(str(self.root))
        self.r.mouse(dict(mouse_id='TEST001',cage='TEST',sex='F',group='KA',weight_g=25))
    def tearDown(self):
        self.r.request_stop('test cleanup')
        if self.r.worker:self.r.worker.join(12)
        self.r.close_services()
        self.temp.cleanup()
    def payload(self,seconds=2):
        return dict(phase='baseline',mouse_id='TEST001',camera='synthetic',size='640x480',fps=25,
                    test=True,baseline_seconds=seconds,post_seconds=seconds)
    def wait(self,limit=15):
        until=time.monotonic()+limit
        while self.r.active and time.monotonic()<until:time.sleep(.1)
        self.assertIsNone(self.r.active)
    def test_record_management(self):
        self.assertEqual(self.r.capture_settings({'camera':'test'})['size'],'1280x720')
        sid=self.r.start(self.payload())
        with self.assertRaises(ValueError):
            self.r.manage_session(dict(session_id=sid,action='delete',reason='running'))
        self.wait()
        with self.assertRaises(ValueError):
            self.r.manage_session(dict(session_id=sid,action='edit',reason='waiting'))
        self.r.event(dict(session_id=sid,kind='abandon',note='management test'))
        self.r.session(sid)['test']=False  # Exercise experiment export filtering.
        folder=Path(self.r.session(sid)['folder'])
        original=(folder/'baseline.mkv').read_bytes()
        with self.assertRaises(ValueError):
            self.r.manage_session(dict(session_id=sid,action='edit',mouse={'weight_g':-1},reason='bad'))
        self.assertEqual(self.r.session(sid)['mouse']['weight_g'],25)
        self.r.manage_session(dict(session_id=sid,action='edit',mouse={'weight_g':20,'notes':'corrected'},reason='weighing correction'))
        current=self.r.session(sid)
        self.assertEqual(current['mouse']['volume_ml'],.25)
        self.assertEqual(self.r.db['mice'][0]['weight_g'],25)
        self.assertEqual(current['revisions'][0]['before']['mouse']['weight_g'],25)
        self.r.manage_session(dict(session_id=sid,action='delete',reason='test record'))
        self.assertTrue(self.r.session(sid)['deleted_at'])
        with zipfile.ZipFile(io.BytesIO(export_workbook(self.r.state()))) as z:
            self.assertNotIn(sid.encode(),z.read('xl/worksheets/sheet2.xml'))
        with self.assertRaises(ValueError):
            self.r.manage_session(dict(session_id=sid,action='edit',reason='deleted'))
        self.r.close_services()
        restored=Recorder(self.root/'state');self.addCleanup(restored.close_services)
        self.assertTrue(restored.session(sid)['deleted_at'])
        restored.manage_session(dict(session_id=sid,action='restore',reason='undo deletion'))
        self.assertNotIn('deleted_at',restored.session(sid))
        with zipfile.ZipFile(io.BytesIO(export_workbook(restored.state()))) as z:
            self.assertIn(sid.encode(),z.read('xl/worksheets/sheet2.xml'))
        self.assertEqual((folder/'baseline.mkv').read_bytes(),original)
        self.assertEqual(len(restored.session(sid)['revisions']),3)

    def test_volume(self):
        for w,v in [(20,.25),(25,.3125),(30,.375)]:self.assertAlmostEqual(volume(w),v)
        for w in [0,-1,float('nan'),float('inf'),'']:
            with self.assertRaises(ValueError):volume(w)
    def test_two_phases_and_persistence(self):
        sid=self.r.start(self.payload())
        with self.assertRaises(ValueError):self.r.start(self.payload())
        self.wait()
        s=self.r.session(sid)
        self.assertEqual(s['status'],'awaiting_injection',s.get('error'))
        with self.assertRaises(ValueError):self.r.start(dict(phase='post',session_id=sid))
        with self.assertRaises(ValueError):self.r.event(dict(kind='injection',session_id=sid,actual_volume_ml=0))
        self.r.event(dict(kind='injection',session_id=sid,actual_volume_ml=.3125))
        with self.assertRaises(ValueError):self.r.event(dict(kind='injection',session_id=sid,actual_volume_ml=.3125))
        self.r.start(dict(phase='post',session_id=sid));self.wait()
        self.assertEqual(s['status'],'completed',s.get('error'))
        self.assertEqual(len(list(Path(s['folder']).glob('*.mkv'))),2)
        self.assertEqual(len(list(Path(s['folder']).glob('*.mp4'))),0)
        self.assertTrue(all(q['mp4_status']=='queued' for q in s['phases'].values()))
        self.r.queue_action({'action':'process'})
        until=time.monotonic()+20
        while any(q['mp4_status']!='verified' for q in s['phases'].values()) and time.monotonic()<until:time.sleep(.1)
        self.assertEqual(len(list(Path(s['folder']).glob('*.mp4'))),2)
        for phase in s['phases'].values():
            self.assertEqual(phase['mp4_status'],'verified',phase)
            self.assertEqual(phase['qc'],'TIMING_CHECKS_PASSED',phase)
        self.assertGreater(s['phases']['post']['frames'],30)
        self.r.close_services()
        restored=Recorder(self.root/'state');self.addCleanup(restored.close_services)
        self.assertEqual(restored.session(sid)['status'],'completed')
    def test_interruption_and_recovery(self):
        sid=self.r.start(self.payload(30));time.sleep(1.8)
        self.r.request_stop('operator test stop');self.wait()
        self.assertEqual(self.r.session(sid)['status'],'interrupted')
        self.r.session(sid)['status']='recording';self.r.persist()
        self.r.close_services()
        recovered=Recorder(self.root/'state');self.addCleanup(recovered.close_services)
        self.assertEqual(recovered.session(sid)['status'],'interrupted')
    def test_encoder_failure_is_not_completion(self):
        sid=self.r.start(self.payload(30));time.sleep(1.2)
        self.r.proc.kill();self.wait()
        self.assertEqual(self.r.session(sid)['status'],'interrupted')
        self.assertNotEqual(self.r.session(sid)['phases']['baseline']['returncode'],0)
    def test_disk_preflight_and_missing_mouse(self):
        with patch('app.shutil.disk_usage',return_value=SimpleNamespace(free=0)):
            with self.assertRaises(ValueError):self.r.start(self.payload())
        self.assertEqual(self.r.db['sessions'],[])
        bad=self.payload();bad['mouse_id']='UNREGISTERED'
        with self.assertRaises(StopIteration):self.r.start(bad)
        self.assertIsNone(self.r.active)
    def test_export_and_literal_id(self):
        self.r.mouse(dict(mouse_id='=1+1',cage='A',sex='M',group='Saline',weight_g=20))
        b=export_workbook(self.r.state())
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            self.assertIsNone(z.testzip())
            r=E.fromstring(z.read('xl/worksheets/sheet1.xml'))
            ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            c=r.find('.//s:c[@r="A9"]',ns)
            self.assertEqual(c.attrib['t'],'inlineStr')
            self.assertEqual(c.find('s:is/s:t',ns).text,'=1+1')

if __name__=='__main__':unittest.main(verbosity=2)
