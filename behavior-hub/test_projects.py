import unittest,tempfile,time,json,shutil,threading
from pathlib import Path
from unittest.mock import patch
from app import Recorder,export_log
from media_quality import ProcessingCancelled

class ProjectTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name);self.r=Recorder(self.root/'state',True)
  self.r.project_action(dict(action='new',name='Mixed behaviors',folder=str(self.root)))
 def tearDown(self):
  self.r.request_stop('test cleanup')
  if self.r.worker:self.r.worker.join(12)
  self.r.close_services();self.temp.cleanup()
 def payload(self,assay='OFT',seconds=1):return dict(camera='synthetic',size='160x120',fps=25,assay=assay,duration_seconds=seconds,test=True)
 def wait(self,predicate,seconds=20):
  end=time.monotonic()+seconds
  while not predicate() and time.monotonic()<end:time.sleep(.05)
  self.assertTrue(predicate())
 def record(self,assay='OFT'):
  sid=self.r.start(self.payload(assay));self.wait(lambda:not self.r.active);return sid
 def test_consecutive_recordings_then_idle_conversion(self):
  sid=self.record();s=self.r.session(sid);self.assertEqual(s['status'],'completed')
  self.assertEqual(s['phases']['recording']['mp4_status'],'queued');self.assertFalse(list(Path(s['folder']).glob('*.mp4')))
  second=self.record('NPR');self.assertNotEqual(sid,second);self.assertRegex(sid,r'^\d{8}_\d{6}_OFT_001_[0-9a-f]{8}$')
  self.r.queue_action({'action':'process'});self.wait(lambda:all(x['phases']['recording'].get('mp4_status')=='verified' for x in self.r.db['sessions']))
  self.assertEqual(s['phases']['recording']['decoded_frames'],25)
 def test_moved_project_relative_paths_and_settings(self):
  sid=self.record();old=self.r.project_file;copy=self.root/'moved';shutil.copytree(old.parent,copy,ignore=shutil.ignore_patterns('*.lock'))
  raw=json.loads((copy/old.name).read_text());self.assertFalse(Path(raw['sessions'][0]['folder']).is_absolute());self.assertEqual(raw['output'],'data')
  self.r.project_action(dict(action='open',path=str(copy/old.name)));s=self.r.session(sid)
  self.assertTrue(Path(s['folder']).is_relative_to(copy));self.assertEqual(self.r.db['settings']['assay'],'OFT')
  self.r.queue_action({'action':'process'});self.wait(lambda:s['phases']['recording'].get('mp4_status')=='verified')
  self.assertFalse(list(old.parent.rglob('*.mp4')))
 def test_queue_preempted_before_recording(self):
  sid=self.record();entered=threading.Event()
  def slow(*args,**kwargs):
   cancel=args[-1];entered.set()
   while not cancel.wait(.02):pass
   raise ProcessingCancelled()
  with patch('background_jobs.finalize_video',side_effect=slow):
   self.r.queue_action({'action':'process'});self.assertTrue(entered.wait(5));before=time.monotonic()
   second=self.r.start(self.payload('FST'));self.assertLess(time.monotonic()-before,3)
   self.assertIsNone(self.r.processing_id);self.assertEqual(self.r.session(sid)['phases']['recording']['mp4_status'],'queued')
   self.wait(lambda:not self.r.active)
 def test_project_isolation_lock_and_path_escape(self):
  sid=self.record();old=self.r.project_file
  other=Recorder(self.root/'other_state',True)
  try:
   with self.assertRaises(OSError):other.project_action({'action':'open','path':str(old)})
  finally:other.close_services()
  self.r.project_action({'action':'new','folder':str(self.root),'name':'second'})
  self.assertEqual(self.r.db['sessions'],[])
  self.r.project_action({'action':'open','path':str(old)});self.assertEqual(self.r.db['sessions'][0]['id'],sid)
  bad=self.root/'bad.project.json';data=json.loads(old.read_text());data['sessions'][0]['folder']='../outside';bad.write_text(json.dumps(data))
  with self.assertRaises(ValueError):self.r.project_action({'action':'open','path':str(bad)})
  self.assertEqual(self.r.project_file,old)
 def test_failure_retry_and_restart_queue(self):
  sid=self.record();s=self.r.session(sid);path=Path(s['folder'])/s['phases']['recording']['file'];original=path.read_bytes();path.write_bytes(b'bad')
  self.r.queue_action({'action':'process'});self.wait(lambda:s['phases']['recording']['mp4_status']=='failed')
  self.assertEqual(s['status'],'completed');path.write_bytes(original)
  self.r.queue_action({'action':'retry'});self.assertEqual(s['phases']['recording']['mp4_status'],'queued')
  self.r.queue_action({'action':'pause'});file=self.r.project_file;self.r.close_services()
  self.r=Recorder(self.root/'state',True);self.r.project_action({'action':'open','path':str(file)})
  self.assertEqual(self.r.session(sid)['phases']['recording']['mp4_status'],'queued')
  self.r.queue_action({'action':'process'});self.wait(lambda:self.r.session(sid)['phases']['recording']['mp4_status']=='verified')
 def test_real_subprocess_cancel_and_saved_settings(self):
  from media_quality import run_cancellable
  from app import ffmpeg,FLAGS
  cancel=threading.Event();timer=threading.Timer(.3,cancel.set);timer.start();before=time.monotonic()
  try:
   with self.assertRaises(ProcessingCancelled):run_cancellable([ffmpeg(),'-re','-f','lavfi','-i','testsrc2=size=160x120:rate=25','-t','30','-f','null','-'],FLAGS,40,cancel)
   self.assertLess(time.monotonic()-before,3)
  finally:timer.cancel()
  settings={'assay':'FST','duration_seconds':300,'test':False,'capture':self.payload()}
  self.r.project_action({'action':'save','settings':settings});file=self.r.project_file
  self.r.project_action({'action':'legacy'});self.r.project_action({'action':'open','path':str(file)})
  self.assertEqual(self.r.db['settings'],settings)
 def test_automatic_idle_and_pause(self):
  sid=self.record();self.r.queue_action({'action':'pause'});self.r.last_activity=0;time.sleep(.5)
  self.assertEqual(self.r.session(sid)['phases']['recording']['mp4_status'],'queued')
  self.r.queue_action({'action':'resume'});self.r.last_activity=0
  self.wait(lambda:self.r.session(sid)['phases']['recording']['mp4_status']=='verified')
if __name__=='__main__':unittest.main(verbosity=2)
