"""Focus safety/persistence/transport tests; no physical camera required."""
import copy
import json
from pathlib import Path
import subprocess
import tempfile
import time
import unittest
from unittest.mock import patch

from app import Recorder, ROOT, ffmpeg, FLAGS, device_details
from focus_capture import focus_request, validate_report, FocusBridge


def report(request=None):
    return dict(protocol=1, source='active_capture_filter', minimum=0, maximum=250,
                step=5, default=0, capabilities=3, value=40, flags=2,
                verified=True, camera_id='@device_test', requested=request)


class FocusTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.r = Recorder(self.root/'state', synthetic=True)
        self.r.project_action(dict(action='new',name='Focus test',folder=str(self.root)))
        self.c = dict(camera='C920',camera_id='@device_test',size='160x120',fps=25,
                      input_format='mjpeg',focus={'mode':'manual','value':40})

    def tearDown(self):
        self.r.request_stop('test cleanup')
        if self.r.worker: self.r.worker.join(15)
        self.r.close_services()
        self.temp.cleanup()

    def test_invalid_manual_and_readback_fail_closed(self):
        for v in [True, None, '40', float('nan'), float('inf'), 2.5]:
            with self.assertRaises(ValueError): focus_request(dict(mode='manual',value=v))
        request=focus_request(self.c['focus'])
        validate_report(report(),request)
        for key,value in [('verified',False),('source','other_camera_instance'),('flags',1),
                          ('flags',3),('value',45),('step',0),('maximum',30),('capabilities',1)]:
            bad=report();bad[key]=value
            with self.assertRaises(ValueError): validate_report(bad,request)
        with self.assertRaises(ValueError): validate_report(report(),dict(mode='manual',value=41))

    def test_failed_focus_never_launches_encoder_or_creates_video(self):
        self.r.synthetic=False
        for message in ['unsupported focus','readback mismatch','camera unplugged','verification timeout']:
            with patch.object(self.r,'focus_bridge_factory',side_effect=ValueError(message)), patch('app.subprocess.Popen') as encoder:
                with self.assertRaisesRegex(ValueError,'Recording blocked'):
                    self.r.start(dict(self.c,assay='OFT',duration_seconds=2,test=True))
                encoder.assert_not_called()
            self.assertIsNone(self.r.active)
            s=self.r.db['sessions'][-1]
            self.assertEqual(s['status'],'failed')
            self.assertFalse(s['focus']['verified'])
            self.assertFalse(list(Path(s['folder']).glob('*.mkv')))

    def test_saved_profile_is_device_scoped_and_requires_live_verification(self):
        self.r.synthetic=False
        with self.assertRaises(ValueError): self.r.focus_action(dict(self.c,action='save'))
        self.r.active='preview';self.r.live_capture=copy.deepcopy(self.c);self.r.focus_status=report(self.c['focus'])
        self.r.focus_action(dict(self.c,action='save'))
        # Closing the camera makes its previous verification insufficient to save.
        self.r.active=None;self.r.live_capture=None
        with self.assertRaises(ValueError): self.r.focus_action(dict(self.c,action='save'))
        data=dict(self.c);data.pop('focus')
        self.assertEqual(self.r.capture_settings(data)['focus'],self.c['focus'])
        data['camera_id']='@device_other'
        self.assertEqual(self.r.capture_settings(data)['focus']['mode'],'device')
        saved=json.loads(self.r.project_file.read_text())
        self.assertEqual(saved['focus_profiles']['@device_test']['focus'],self.c['focus'])
        path=str(self.r.project_file)
        self.r.project_action({'action':'legacy'})
        self.r.project_action({'action':'open','path':path})
        self.assertEqual(self.r.db['focus_profiles'],saved['focus_profiles'])

    def test_focus_controls_locked_during_recording_and_simulation_not_verified(self):
        self.r.active='recording-id'
        for kind in ['read','apply','save']:
            with self.assertRaisesRegex(ValueError,'locked'): self.r.focus_action(dict(self.c,action=kind))
        self.r.active=None
        with self.assertRaisesRegex(ValueError,'Simulated'): self.r.focus_action(dict(self.c,action='read'))

    def test_duplicate_names_keep_unique_device_ids(self):
        stderr=b'"C920" (video)\n  Alternative name "@device_A"\n"C920" (video)\n  Alternative name "@device_B"\n"Mic" (audio)\n  Alternative name "@device_M"\n'
        with patch('app.run_ff',return_value=subprocess.CompletedProcess([],0,stderr=stderr)):
            self.assertEqual(device_details(),[dict(name='C920',id='@device_A'),dict(name='C920',id='@device_B')])

    def test_nut_relay_preserves_variable_frame_timestamps(self):
        # Exercise the actual relay with synthetic variable-rate NUT. This verifies
        # transport, not the physical camera or native IAMCameraControl calls.
        import socket, threading
        source=self.root/'source.nut'; dest=self.root/'received.nut'
        cmd=[ffmpeg(),'-v','error','-f','lavfi','-i','testsrc2=size=160x120:rate=25',
             '-vf',"select='not(eq(mod(n,7),0))'",'-frames:v','25','-fps_mode','passthrough',
             '-c:v','rawvideo','-f','nut',str(source)]
        subprocess.run(cmd,check=True,capture_output=True,creationflags=FLAGS)
        # A real child process supplies bytes just as the native capture helper does.
        bridge=FocusBridge.__new__(FocusBridge)
        bridge.proc=subprocess.Popen([ffmpeg(),'-v','error','-i',str(source),'-c:v','copy','-f','nut','pipe:1'],
                                     stdout=subprocess.PIPE,stderr=subprocess.PIPE,creationflags=FLAGS,bufsize=0)
        bridge.closed=threading.Event();bridge.listener=None;bridge.connection=None;bridge.relay_thread=None
        bridge.reader=threading.Thread(target=lambda:bridge.proc.stderr.read(),daemon=True);bridge.reader.start()
        try:
            subprocess.run([ffmpeg(),'-v','error',*bridge.input_args(),'-c:v','copy',str(dest)],
                           check=True,capture_output=True,timeout=15,creationflags=FLAGS)
        finally: bridge.close()
        def hashes(path):
            p=subprocess.run([ffmpeg(),'-v','error','-i',str(path),'-f','framemd5','-'],capture_output=True,check=True,creationflags=FLAGS)
            return [x for x in p.stdout.decode().splitlines() if not x.startswith('#')]
        self.assertEqual(hashes(source),hashes(dest))


if __name__=='__main__': unittest.main(verbosity=2)
