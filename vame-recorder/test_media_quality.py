import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from app import ffmpeg, run_ff
from media_quality import finalize_video,timing_metrics

class MediaQualityTests(unittest.TestCase):
    def test_gap_and_deficit(self):
        q=timing_metrics([0,.04,.08,.20],25,.24)
        self.assertEqual(q['qc'],'REVIEW')
        self.assertEqual(q['estimated_missing_from_gaps'],2)
        self.assertEqual(q['frame_count_deficit'],2)
        self.assertAlmostEqual(q['max_gap_ms'],120)
    def test_regular_and_bad_timestamps(self):
        self.assertEqual(timing_metrics([i/25 for i in range(50)],25,2)['qc'],'TIMING_CHECKS_PASSED')
        self.assertEqual(timing_metrics([0,0,.04],25,.12)['non_increasing_timestamps'],1)
    def test_real_remux_preserves_variable_timestamps(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'baseline.mkv'
            r=run_ff(['-f','lavfi','-i','testsrc2=size=160x120:rate=25','-t','2','-vf',"select=not(eq(n\\,10))",'-fps_mode','passthrough','-c:v','libx264',str(p)])
            self.assertEqual(r.returncode,0,r.stderr)
            original=p.read_bytes()
            q=finalize_video(ffmpeg(),p,25,2)
            self.assertEqual(q['mp4_status'],'verified',q)
            self.assertEqual(q['decoded_frames'],49)
            self.assertEqual(q['gaps_over_threshold'],1)
            self.assertEqual(q['qc'],'REVIEW')
            self.assertEqual(original,p.read_bytes())
            self.assertTrue(p.with_suffix('.mp4').exists())
            self.assertTrue(p.with_name('baseline_frame_timing.csv').exists())
            # Existing MP4 must not be overwritten.
            before=p.with_suffix('.mp4').read_bytes()
            again=finalize_video(ffmpeg(),p,25,2)
            self.assertEqual(again['mp4_status'],'verified')
            self.assertEqual(before,p.with_suffix('.mp4').read_bytes())
    def test_decode_failure_retains_source(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'bad.mkv';p.write_bytes(b'broken test video')
            q=finalize_video(ffmpeg(),p,25,2)
            self.assertEqual(q['mp4_status'],'failed')
            self.assertEqual(q['qc'],'REVIEW')
            self.assertEqual(p.read_bytes(),b'broken test video')
    def test_conversion_failure_is_not_verified(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'source.mkv';p.write_bytes(b'original')
            with patch('media_quality.decode',return_value=[(0,'hash'),(.04,'hash')]),patch('media_quality.run_cancellable',side_effect=OSError('disk full')):
                q=finalize_video(ffmpeg(),p,25,.08)
            self.assertEqual(q['mp4_status'],'failed')
            self.assertIn('disk full',q['mp4_error'])
            self.assertTrue(p.exists())

if __name__=='__main__':unittest.main(verbosity=2)
