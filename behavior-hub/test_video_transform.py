import io
import subprocess
import tempfile
import time
import unittest
import json
from pathlib import Path
from PIL import Image
from app import Recorder, ffmpeg, FLAGS
from video_transform import transform_settings, framing_filter

class TransformTests(unittest.TestCase):
    def test_framing_updates_keep_camera_and_frames_alive(self):
        with tempfile.TemporaryDirectory() as temp:
            r=Recorder(Path(temp)/'state',True)
            c=dict(camera='synthetic',size='160x120',fps=25)
            try:
                # Also works with no project open, as in the reported issue.
                r.preview_start(c)
                end=time.monotonic()+15
                while not r.preview and time.monotonic()<end:time.sleep(.05)
                self.assertTrue(r.preview)
                proc=r.proc; worker=r.worker; before=r.preview_at
                for i in range(30):
                    r.preview_start(dict(c,transform=dict(zoom=1+i/10,pan_x=(i%3)-1,flip_h=bool(i%2))))
                    self.assertIs(r.proc,proc);self.assertIs(r.worker,worker)
                    self.assertEqual(r.active,'preview');self.assertTrue(r.preview)
                end=time.monotonic()+3
                while r.preview_at<=before and time.monotonic()<end:time.sleep(.02)
                self.assertGreater(r.preview_at,before)
                self.assertIsNone(proc.poll())
                self.assertEqual(r.live_capture['transform']['zoom'],3.9)
                with self.assertRaises(ValueError):r.preview_start(dict(c,transform={'zoom':8}))
                self.assertIs(r.proc,proc)
                # Actual capture changes must still restart/reconfigure the source.
                r.preview_start(dict(c,size='320x240'))
                self.assertIsNot(r.proc,proc)
            finally:
                r.request_stop('cleanup')
                if r.worker:r.worker.join(15)
                r.close_services()

    def test_browser_geometry_matches_encoder(self):
        # Execute the real browser helper, compare fractional zoom/edge crops and
        # all flip combinations against FFmpeg's recorded crop coordinates.
        cases=[dict(zoom=z,pan_x=x,pan_y=-x,flip_h=h,flip_v=v)
               for z in [1,1.3,2.7,4] for x in [-1,-.3,0,.7,1]
               for h in [False,True] for v in [False,True]]
        script="""const fs=require('fs');const h=fs.readFileSync('index.html','utf8');
        eval(h.slice(h.indexOf('function previewFraming('),h.indexOf('function drawPreviewFraming(')));
        const cases=JSON.parse(fs.readFileSync(0,'utf8'));
        console.log(JSON.stringify(cases.map(t=>previewFraming('1280x720',t))));"""
        result=subprocess.run(['node','-e',script],input=json.dumps(cases),text=True,capture_output=True,check=True)
        for t,f in zip(cases,json.loads(result.stdout)):
            vf=framing_filter('1280x720',t)
            if t['zoom']>1:self.assertTrue(vf.startswith(f"crop={f['w']}:{f['h']}:{f['x']}:{f['y']}"))
            for px,py in [(f['x'],f['y']),(f['x']+f['w'],f['y']+f['h'])]:
                out_x=f['sx']*px/1280+f['tx'];out_y=f['sy']*py/720+f['ty']
                edge=int(px!=f['x'])
                self.assertAlmostEqual(out_x,1-edge if t['flip_h'] else edge)
                self.assertAlmostEqual(out_y,1-edge if t['flip_v'] else edge)

    def test_validation(self):
        for value in [[],{'zoom':0},{'zoom':5},{'zoom':True},{'pan_x':float('nan')},{'pan_y':2},{'flip_h':1}]:
            with self.assertRaises(ValueError):transform_settings(value)
        self.assertEqual(framing_filter('128x96'), 'null')
        self.assertEqual(transform_settings({'zoom':1,'pan_x':1})['pan_x'],0)

    def test_real_preview_recording_and_restore(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'quadrants.ppm'
            colors=[(255,0,0),(0,255,0),(0,0,255),(255,255,0)]
            source.write_bytes(b'P6\n128 96\n255\n'+b''.join(bytes(colors[(y>=48)*2+(x>=64)]) for y in range(96) for x in range(128)))
            r=Recorder(root/'state',True)
            r.project_action(dict(action='new',name='Transform test',folder=str(root)))
            r.input_args=lambda c:['-re','-loop','1','-framerate','25','-i',str(source)]
            def wait(predicate):
                end=time.monotonic()+20
                while not predicate() and time.monotonic()<end:time.sleep(.05)
                self.assertTrue(predicate())
            def color_near(image,point,color):
                value=image.convert('RGB').getpixel(point)
                self.assertTrue(all(abs(a-b)<35 for a,b in zip(value,color)),(value,color))
            try:
                for transform,expected in [({'flip_h':True},colors[1]),({'flip_v':True},colors[2]),({'flip_h':True,'flip_v':True},colors[3]),({'zoom':2,'pan_x':1,'pan_y':1},colors[3])]:
                    c=dict(camera='synthetic',size='128x96',fps=25,transform=transform)
                    r.preview_start(c);wait(lambda:bool(r.preview))
                    preview=Image.open(io.BytesIO(r.preview));self.assertEqual(preview.size,(128,96));color_near(preview,(16,12),colors[0])
                    sid=r.start(dict(c,assay='OFT',test=True,duration_seconds=1));wait(lambda:not r.active)
                    s=r.session(sid);video=Path(s['folder'])/s['phases']['recording']['file']
                    result=subprocess.run([ffmpeg(),'-v','error','-i',str(video),'-frames:v','1','-f','image2pipe','-vcodec','png','pipe:1'],capture_output=True,check=True,creationflags=FLAGS)
                    frame=Image.open(io.BytesIO(result.stdout));self.assertEqual(frame.size,(128,96));color_near(frame,(16,12),expected)
                    r.queue_action({'action':'process'});wait(lambda:s['phases']['recording'].get('mp4_status')=='verified');wait(lambda:not r.processing_id)
                    self.assertEqual(s['capture']['transform'],transform_settings(transform))
                file=r.project_file;r.close_services();r=Recorder(root/'state',True)
                r.project_action(dict(action='open',path=str(file)))
                self.assertEqual(r.db['settings']['capture']['transform']['zoom'],2)
            finally:
                r.request_stop('cleanup')
                if r.worker:r.worker.join(15)
                r.close_services()

if __name__=='__main__':unittest.main(verbosity=2)
