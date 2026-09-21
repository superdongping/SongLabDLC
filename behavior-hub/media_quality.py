"""Post-capture timing QC and lossless container conversion. Originals are never deleted."""
import csv
from fractions import Fraction
import math
from pathlib import Path
import subprocess


def timing_metrics(pts, fps, target):
    gaps = [b-a for a,b in zip(pts, pts[1:])]
    period = 1 / fps
    expected = round(target * fps)
    span = pts[-1]-pts[0] if len(pts)>1 else 0
    long_gaps = [g for g in gaps if g > period * 1.5 + 1e-6]
    bad_order = sum(g <= 0 for g in gaps)
    reasons = []
    if abs(len(pts)-expected) / max(expected, 1) > .005:
        reasons.append('Frame count differs from the target by more than 0.5%.')
    if long_gaps: reasons.append('Frame intervals exceed 1.5 times the target interval.')
    if bad_order: reasons.append('Duplicate or non-increasing frame timestamps.')
    if len(pts)<2: reasons.append('Too few frames to measure timing.')
    return dict(decoded_frames=len(pts), expected_frames=expected,
                observed_fps=(len(pts)-1)/span if span>0 else None,
                max_gap_ms=max(gaps)*1000 if gaps else None,
                gaps_over_threshold=len(long_gaps), gap_threshold_ms=period*1500,
                estimated_missing_from_gaps=sum(max(0,math.floor(g/period+.5)-1) for g in long_gaps),
                non_increasing_timestamps=bad_order, frame_count_deficit=max(0,expected-len(pts)),
                qc='REVIEW' if reasons else 'TIMING_CHECKS_PASSED', qc_reasons=reasons)


def decode(ffmpeg, path, flags, timeout):
    result = subprocess.run([ffmpeg,'-hide_banner','-v','error','-xerror','-i',str(path),
                             '-map','0:v:0','-an','-fps_mode','passthrough','-enc_time_base','1:1000000','-f','framemd5','-'],
                            capture_output=True,timeout=timeout,creationflags=flags)
    if result.returncode:
        raise RuntimeError('Decode failed: '+result.stderr.decode('utf-8','replace')[-1200:])
    tb = None
    rows = []
    for line in result.stdout.decode('ascii').splitlines():
        if line.startswith('#tb 0:'): tb=Fraction(line.split(':',1)[1].strip())
        elif line and not line.startswith('#'):
            parts=[v.strip() for v in line.split(',')]
            rows.append((float(int(parts[2])*tb),parts[5]))
    if not rows: raise RuntimeError('No decoded frames.')
    return rows


def finalize_video(ffmpeg, source, fps, target, flags=0):
    result = dict(qc='REVIEW',qc_reasons=[],mp4_status='failed',mp4_file=None)
    timeout = max(120, target*2)
    try:
        original=decode(ffmpeg,source,flags,timeout)
        result.update(timing_metrics([x[0] for x in original],fps,target))
        csv_path=source.with_name(source.stem+'_frame_timing.csv')
        with csv_path.open('w',newline='',encoding='utf-8') as out:
            writer=csv.writer(out)
            writer.writerow(['frame_index','relative_pts_seconds','interval_seconds','decoded_frame_md5'])
            for i,(pts,checksum) in enumerate(original):
                writer.writerow([i,format(pts-original[0][0],'.9f'),format(pts-original[i-1][0],'.9f') if i else '',checksum])
        result['timing_file']=csv_path.name
        final=source.with_suffix('.mp4')
        partial=source.with_name(source.stem+'.pending.mp4')
        remux=subprocess.run([ffmpeg,'-hide_banner','-v','error','-i',str(source),'-map','0:v:0',
                              '-an','-c:v','copy','-movflags','+faststart','-n',str(partial)],
                             capture_output=True,timeout=timeout,creationflags=flags)
        if remux.returncode: raise RuntimeError('MP4 conversion failed: '+remux.stderr.decode('utf-8','replace')[-1200:])
        converted=decode(ffmpeg,partial,flags,timeout)
        if len(original)!=len(converted) or any(a[1]!=b[1] for a,b in zip(original,converted)):
            raise RuntimeError('MP4 decoded frame count or frame content differs from MKV.')
        error=max(abs((a[0]-original[0][0])-(b[0]-converted[0][0])) for a,b in zip(original,converted))
        if error>.0011: raise RuntimeError('MP4 relative timestamps differ by more than 1.1 ms.')
        if final.exists(): raise RuntimeError('MP4 destination already exists; it was not overwritten.')
        partial.rename(final)
        result.update(mp4_status='verified',mp4_file=final.name,mp4_bytes=final.stat().st_size,
                      mp4_max_relative_timestamp_error_ms=error*1000,
                      mp4_verification='All decoded frame hashes and counts match MKV; relative PTS within 1.1 ms.')
    except Exception as exc:
        result['mp4_error']=str(exc)
        result['qc']='REVIEW'
        result['qc_reasons'].append('Post-capture validation/conversion failed. MKV retained; pending MP4 is not verified.')
    return result
