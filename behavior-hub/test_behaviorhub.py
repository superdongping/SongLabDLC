import unittest,tempfile,time,json,io,csv
from pathlib import Path
from app import Recorder,export_log
class BehaviorHubTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name);self.r=Recorder(self.root/'state',True);self.r.project_action({'action':'new','name':'Test project','folder':str(self.root)})
 def tearDown(self):
  self.r.request_stop('cleanup')
  if self.r.worker:self.r.worker.join(30)
  self.r.close_services();self.temp.cleanup()
 def payload(self,assay='OFT',seconds=1):return dict(camera='synthetic',size='160x120',fps=25,assay=assay,duration_seconds=seconds,test=True)
 def wait(self):
  end=time.monotonic()+30
  while self.r.active and time.monotonic()<end:time.sleep(.1)
  self.assertIsNone(self.r.active)
 def test_six_assays_optional_fields_and_files(self):
  for assay in ['OFT','NPR','ZERO_MAZE','Y_MAZE','FST','TST']:
   sid=self.r.start(self.payload(assay));self.wait();s=self.r.session(sid)
   self.assertEqual(s['status'],'completed',s)
   self.assertEqual(s['mouse']['mouse_id'],'');self.assertIsNone(s['mouse']['weight_g'])
   q=s['phases']['recording'];self.assertEqual(q['mp4_status'],'queued',q)
   self.r.queue_action({'action':'process'})
   end=time.monotonic()+15
   while q['mp4_status'] in ('queued','processing') and time.monotonic()<end:time.sleep(.1)
   self.assertEqual(q['mp4_status'],'verified',q)
   self.assertEqual(q['qc'],'TIMING_CHECKS_PASSED',q)
   self.assertTrue((Path(s['folder'])/(sid+'.mkv')).exists());self.assertTrue((self.r.project_file.parent/q['mp4_relative']).exists())
   self.assertEqual(Path(s['folder']).parent.name,assay)
  self.assertEqual(len({s['id'] for s in self.r.db['sessions']}),6)
 def test_management_and_restart(self):
  sid=self.r.start(self.payload());self.wait();s=self.r.session(sid);video=Path(s['folder'])/(sid+'.mkv');before=video.read_bytes()
  self.r.manage_session(dict(session_id=sid,action='edit',mouse={'mouse_id':'=formula','notes':'correction'},reason='test'))
  self.assertIn("'=formula",export_log(self.r.state()).decode('utf-8-sig'))
  self.r.manage_session(dict(session_id=sid,action='delete',reason='hide'))
  self.assertNotIn(sid,export_log(self.r.state()).decode('utf-8-sig'))
  project=self.r.project_file;self.r.close_services();r=Recorder(self.root/'state');r.project_action({'action':'open','path':str(project)});self.r=r;r.manage_session(dict(session_id=sid,action='restore',reason='restore'))
  self.assertIn(sid,export_log(r.state()).decode('utf-8-sig'));self.assertEqual(before,video.read_bytes())
 def test_interrupted_and_protected(self):
  sid=self.r.start(self.payload(seconds=30));time.sleep(1.5)
  with self.assertRaises(ValueError):self.r.start(self.payload())
  with self.assertRaises(ValueError):self.r.manage_session(dict(session_id=sid,action='delete',reason='active'))
  self.r.request_stop('test early stop');self.wait();self.assertEqual(self.r.session(sid)['status'],'interrupted')
 def test_presets_and_validation(self):
  p=json.loads(Path('presets.json').read_text());self.assertEqual(p['seconds'],dict(OFT=360,NPR=360,ZERO_MAZE=360,Y_MAZE=480,FST=300,TST=360))
  self.assertEqual(self.r.capture_settings({'camera':'test'})['size'],'1280x720')
  for key,value in [('weight_g',-1),('duration_seconds',0),('assay','invalid')]:
   x=self.payload();x[key]=value
   with self.assertRaises(ValueError):self.r.start(x)
if __name__=='__main__':unittest.main(verbosity=2)
