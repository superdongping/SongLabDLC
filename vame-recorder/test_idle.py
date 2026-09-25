import unittest,time,threading,tempfile
from pathlib import Path
from unittest.mock import patch
from test_recorder import RecorderTests
from app import Recorder,ffmpeg,FLAGS
from media_quality import ProcessingCancelled,run_cancellable
class IdleTests(unittest.TestCase):
 setUp=RecorderTests.setUp
 tearDown=RecorderTests.tearDown
 payload=RecorderTests.payload
 wait=RecorderTests.wait
 def until(self,fn):
  end=time.monotonic()+20
  while not fn() and time.monotonic()<end:time.sleep(.05)
  self.assertTrue(fn())
 def complete(self):
  sid=self.r.start(self.payload(1));self.wait();self.r.event(dict(session_id=sid,kind='injection',actual_volume_ml=.3125))
  self.r.start(dict(phase='post',session_id=sid));self.wait();return sid
 def test_waiting_injection_and_next_animal_do_not_wait_for_mp4(self):
  sid=self.r.start(self.payload(1));self.wait();self.r.last_activity=0;time.sleep(.7)
  self.assertEqual(self.r.session(sid)['phases']['baseline']['mp4_status'],'queued')
  self.r.event(dict(session_id=sid,kind='injection',actual_volume_ml=.3125));self.r.start(dict(phase='post',session_id=sid));self.wait()
  entered=threading.Event()
  def slow(*args):
   entered.set();args[-1].wait(10);raise ProcessingCancelled()
  with patch('background_jobs.finalize_video',side_effect=slow):
   self.r.queue_action({'action':'process'});self.assertTrue(entered.wait(5));start=time.monotonic()
   nextsid=self.r.start(self.payload(1));self.assertLess(time.monotonic()-start,3);self.wait()
  self.assertNotEqual(sid,nextsid);self.assertEqual(self.r.session(sid)['phases']['baseline']['mp4_status'],'queued')
 def test_queue_restart_failure_retry_and_automatic_idle(self):
  sid=self.complete();self.r.queue_action({'action':'pause'});s=self.r.session(sid);p=Path(s['folder'])/'baseline.mkv';original=p.read_bytes();p.write_bytes(b'bad')
  self.r.close_services();self.r=Recorder(self.root/'state',True);s=self.r.session(sid)
  self.assertEqual(len(self.r.queue_state()['jobs']),2);self.assertTrue(self.r.db['processing']['paused'])
  self.r.queue_action({'action':'process'});self.until(lambda:s['phases']['baseline']['mp4_status']=='failed')
  self.r.queue_action({'action':'pause'});p.write_bytes(original);self.r.queue_action({'action':'retry'});self.r.queue_action({'action':'resume'});self.r.last_activity=0
  self.until(lambda:all(q['mp4_status']=='verified' for q in s['phases'].values()));self.assertEqual(p.read_bytes(),original)
 def test_preview_delivery_and_recording_preview_limit(self):
  self.r.preview_start(self.payload());self.assertEqual(self.r.state()['preview_target_fps'],25)
  self.until(lambda:bool(self.r.preview));seen=set();end=time.monotonic()+2
  while time.monotonic()<end:seen.add(self.r.preview_at);time.sleep(.008)
  self.assertGreater(len(seen),25,'Synthetic preview did not exceed the old low FPS')
  sid=self.r.start(self.payload(1));self.assertEqual(self.r.state()['preview_target_fps'],8);self.wait()
 def test_real_worker_process_cancellation(self):
  cancel=threading.Event();timer=threading.Timer(.3,cancel.set);timer.start();start=time.monotonic()
  try:
   with self.assertRaises(ProcessingCancelled):run_cancellable([ffmpeg(),'-re','-f','lavfi','-i','testsrc2=size=160x120:rate=25','-t','30','-f','null','-'],FLAGS,40,cancel)
   self.assertLess(time.monotonic()-start,3)
  finally:timer.cancel()
del RecorderTests
if __name__=='__main__':unittest.main(verbosity=2)
