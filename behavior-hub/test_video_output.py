import json
import shutil
import tempfile
import time
import unittest
from unittest.mock import patch
from pathlib import Path
from app import Recorder, export_log
from video_output import reserve_mp4, mp4_destination, naming_mode


class VideoOutputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.r = Recorder(self.root/'state', True)
        self.r.project_action(dict(action='new', name='Output test', folder=str(self.root)))

    def tearDown(self):
        self.r.request_stop('cleanup')
        if self.r.worker:
            self.r.worker.join(15)
        self.r.close_services()
        self.temp.cleanup()

    def wait(self, predicate):
        end = time.monotonic()+25
        while not predicate() and time.monotonic() < end:
            time.sleep(.05)
        self.assertTrue(predicate())

    def record(self, mouse='', mode='datetime_mouse'):
        sid = self.r.start(dict(camera='synthetic', size='160x120', fps=25, assay='OFT',
                                duration_seconds=1, test=True, mouse_id=mouse, video_naming=mode))
        self.wait(lambda: not self.r.active)
        return self.r.session(sid)

    def process(self, s, expected='verified'):
        self.r.queue_action({'action':'process'})
        self.wait(lambda: s['phases']['recording'].get('mp4_status') == expected)
        self.wait(lambda: not self.r.processing_id)

    def test_collisions_sanitization_and_escape(self):
        root=self.r.project_file.parent
        one=reserve_mp4(root,[], 'id','20261001_120000','Mouse:/01','datetime_mouse')
        self.assertEqual(one,'MP4/20261001_120000_Mouse__01.mp4')
        sessions=[{'phases':{'recording':{'mp4_relative':one.upper()}}}]
        two=reserve_mp4(root,sessions,'id','20261001_120000','Mouse:/01','datetime_mouse')
        self.assertTrue(two.endswith('_002.mp4'))
        (root/two).write_bytes(b'keep')
        three=reserve_mp4(root,sessions,'id','20261001_120000','Mouse:/01','datetime_mouse')
        self.assertTrue(three.endswith('_003.mp4'))
        self.assertEqual(reserve_mp4(root,[],'id','20261001_120001','','datetime_mouse'),'MP4/20261001_120001.mp4')
        for relative in ['../out.mp4','MP4/../../out.mp4','data/out.mp4','MP4/a.txt']:
            with self.assertRaises(ValueError): mp4_destination(root,relative)
        with self.assertRaises(ValueError): naming_mode('invalid')

    def test_shared_folder_fixed_names_and_move(self):
        s=self.record('Mouse01'); q=s['phases']['recording']; relative=q['mp4_relative']
        self.assertRegex(relative,r'^MP4/TEST_\d{8}_\d{6}_OFT_Mouse01.mp4$')
        self.r.manage_session(dict(session_id=s['id'],action='edit',mouse={'mouse_id':'Corrected'},reason='correction'))
        s=self.r.session(s['id']); self.assertEqual(s['phases']['recording']['mp4_relative'],relative)
        self.process(s)
        with patch('app.os.startfile') as opened:
            actual=self.r.open_recording_video({'session_id':s['id']})
            opened.assert_called_once_with(str(self.r.project_file.parent/relative))
            self.assertEqual(actual,str(self.r.project_file.parent/relative))
        second=self.record('', 'record_id');self.process(second)
        root=self.r.project_file.parent
        self.assertEqual(len(list((root/'MP4').glob('*.mp4'))),2)
        self.assertFalse(list((root/'data').rglob('*.mp4')))
        self.assertEqual(len(list((root/'data').rglob('*.mkv'))),2)
        self.assertIn(relative,export_log(self.r.db).decode('utf-8-sig'))
        old=self.r.project_file; moved=self.root/'moved'
        shutil.copytree(root,moved,ignore=shutil.ignore_patterns('*.lock'))
        self.r.project_action(dict(action='open',path=str(moved/old.name)))
        self.assertEqual(self.r.db['settings']['video_naming'],'datetime_behavior')
        self.assertTrue(mp4_destination(moved,self.r.session(s['id'])['phases']['recording']['mp4_relative']).is_file())

    def test_conflict_never_overwrites_and_retry_after_restart(self):
        s=self.record('M1');q=s['phases']['recording'];p=self.r.project_file.parent/q['mp4_relative']
        p.write_bytes(b'unrelated existing file')
        with patch('app.os.startfile') as opened:
            with self.assertRaises(ValueError):self.r.open_recording_video({'session_id':s['id']})
            opened.assert_not_called()
        self.process(s,'failed');self.assertEqual(p.read_bytes(),b'unrelated existing file')
        p.rename(p.with_suffix('.conflict'))
        project=self.r.project_file;self.r.close_services()
        self.r=Recorder(self.root/'state',True)
        self.r.project_action(dict(action='open',path=str(project)))
        self.r.queue_action({'action':'retry'});s=self.r.session(s['id']);self.process(s)
        self.assertTrue(p.is_file())
        self.assertEqual(s['phases']['recording']['decoded_frames'],25)
        # Simulate publication succeeding just before the queue state was persisted.
        original=p.read_bytes();s['phases']['recording']['mp4_status']='queued'
        self.process(s);self.assertEqual(p.read_bytes(),original)

    def test_custom_names_and_persistent_automatic_ids(self):
        payload=dict(camera='synthetic',size='160x120',fps=25,assay='CUSTOM',duration_seconds=1,test=True)
        with self.assertRaisesRegex(ValueError,'custom behavioral test'):
            self.r.start(payload)
        self.assertEqual(self.r.db['sessions'],[])
        sid=self.r.start(dict(payload,custom_behavior='Social/interaction'))
        self.wait(lambda:not self.r.active)
        s=self.r.session(sid)
        self.assertEqual(s['behavior_name'],'Social/interaction')
        self.assertTrue(s['phases']['recording']['mp4_relative'].endswith('_Social_interaction_Auto_ID01.mp4'))
        self.assertIn('Social/interaction',export_log(self.r.db).decode('utf-8-sig'))
        self.process(s)
        file=self.r.project_file;self.r.close_services()
        self.r=Recorder(self.root/'state',True)
        self.r.project_action(dict(action='open',path=str(file)))
        self.assertEqual(self.r.db['settings']['custom_behavior'],'Social/interaction')
        named=self.record('Animal9')
        self.assertTrue(named['phases']['recording']['mp4_relative'].endswith('_OFT_Animal9.mp4'))
        second=self.record()
        self.assertTrue(second['phases']['recording']['mp4_relative'].endswith('_OFT_Auto_ID02.mp4'))
        self.process(second)

    def test_old_project_sessions_keep_original_destination(self):
        s=self.record();s['phases']['recording'].pop('mp4_relative')
        self.process(s)
        video=Path(s['folder'])/(s['id']+'.mp4')
        self.assertTrue(video.is_file())
        with patch('app.os.startfile') as opened:
            self.r.open_recording_video({'session_id':s['id']})
            opened.assert_called_once_with(str(video.resolve()))
            video.rename(video.with_suffix('.kept'))
            opened.reset_mock()
            with self.assertRaisesRegex(ValueError,'missing'):
                self.r.open_recording_video({'session_id':s['id']})
            opened.assert_not_called()
            s['phases']['recording']['mp4_file']='../outside.mp4'
            with self.assertRaisesRegex(ValueError,'Invalid'):
                self.r.open_recording_video({'session_id':s['id']})


if __name__ == '__main__':
    unittest.main(verbosity=2)
