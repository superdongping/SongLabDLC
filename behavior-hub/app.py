"""Behavior Hub: loopback-only Windows capture service. No cloud services."""
from media_quality import finalize_video
from focus_capture import FocusBridge, focus_request
from video_output import naming_mode, reserve_mp4, mp4_destination
from video_transform import transform_settings, framing_filter
import csv
import argparse
import copy
import ctypes
import datetime as dt
import io
import json
import math
import os
from pathlib import Path
import re
import secrets
import socket
import shutil
import subprocess
import sys
import threading
import time
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs
from urllib.request import urlopen
import zipfile
import xml.etree.ElementTree as ET

# SongScope remains the internal service/state identifier for upgrade compatibility.
ROOT = Path(getattr(sys, '_MEIPASS', Path(__file__).parent))
FLAGS = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
TOKEN = secrets.token_urlsafe(32)
APP_VERSION = '1.2.6-r2'

def now():
    return dt.datetime.now().astimezone().isoformat(timespec='milliseconds')

def ffmpeg():
    bundled = ROOT / 'ffmpeg.exe'
    if bundled.exists():
        return str(bundled)
    import imageio_ffmpeg
    return imageio_ffmpeg.get_ffmpeg_exe()

def run_ff(args, timeout=20):
    return subprocess.run([ffmpeg(), '-hide_banner', *args], capture_output=True,
                          timeout=timeout, creationflags=FLAGS)

def devices():
    return [d['name'] for d in device_details()]

def device_details():
    r = run_ff(['-list_devices', 'true', '-f', 'dshow', '-i', 'dummy'])
    result = []
    pending = None
    for line in r.stderr.decode('utf-8', 'replace').splitlines():
        name = re.search(r'"([^"\r\n]+)" \(video\)', line)
        if name:
            pending = name.group(1)
        elif 'Alternative name' in line and pending:
            match = re.search(r'Alternative name "([^"]+)"', line)
            if match:
                result.append({'name': pending, 'id': match.group(1)})
            pending = None
        elif '(audio)' in line:
            pending = None
    return result

def atomic_json(path, data):
    tmp = path.with_suffix('.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)

def safe_text(v, limit=200):
    v = str(v).strip()
    if len(v) > limit or any(ord(c) < 32 for c in v):
        raise ValueError('Text is too long or contains control characters.')
    return v

def numeric(v, low, high, label):
    x = float(v)
    if not math.isfinite(x) or not low <= x <= high:
        raise ValueError(f'{label} must be between {low} and {high}.')
    return x

sys.path.insert(0, str(Path(__file__).parent / 'vendor'))
from project_store import ProjectSupport, atomic
from background_jobs import BackgroundSupport

class Recorder(ProjectSupport, BackgroundSupport):
    def __init__(self, data_dir, synthetic=False):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = self.data_dir / 'registry.json'
        self.lock = threading.RLock()
        self.command_lock = threading.Lock()
        self.db = {'mice': [], 'sessions': [], 'output': ''}
        if self.db_path.exists():
            self.db = json.loads(self.db_path.read_text(encoding='utf-8'))
        self.proc = None
        self.worker = None
        self.active = None
        self.preview = b''
        self.preview_at = 0
        self.error = ''
        self.stop_reason = None
        self.synthetic = synthetic
        self.frames = 0
        self.media_seconds = 0
        self.started_mono = None
        self.request_mono = 0
        self.last_progress = 0
        self.log_tail = []
        self.focus_status = None
        self.live_capture = None
        self.focus_bridge_factory = FocusBridge
        self.init_projects()
        self.ffmpeg_path=ffmpeg;self.process_flags=FLAGS
        for s in self.db['sessions']:
            if s['status'] in ('starting', 'recording', 'finalizing'):
                s['status'] = 'interrupted'
                s['error'] = 'The service stopped before recording finished. Check the retained video files.'
                self.save_session(s)
        self.persist()
        self.init_background()

    def persist(self):
        if self.project_file:self.project_persist()
        else:
            atomic_json(self.db_path, self.db)
            self.legacy_db=copy.deepcopy(self.db)

    def save_session(self, s):
        saved=copy.deepcopy(s)
        if self.project_file:saved['folder']='.'
        atomic_json(Path(s['folder']) / 'session.json', saved)

    def session(self, sid):
        return next(s for s in self.db['sessions'] if s['id'] == sid)

    def state(self):
        with self.lock:
            d = copy.deepcopy(self.db)
            d.update(self.project_state())
            d.update(processing=self.queue_state())
            d.update(active=self.active, preview=bool(self.proc and self.active == 'preview'),
                     frames=self.frames, media_seconds=self.media_seconds,
                     elapsed=round(time.monotonic() - self.started_mono, 2) if self.started_mono else 0,
                     error=self.error, preview_age=round(time.monotonic() - self.preview_at, 1) if self.preview_at else None,
                     synthetic=self.synthetic, app_version=APP_VERSION,
                     focus_status=copy.deepcopy(self.focus_status), live_capture=copy.deepcopy(self.live_capture))
            return d

    def set_output(self, path):
        if self.project_file:raise ValueError('Project videos are stored inside the project data folder. Create/open another project to change location.')
        p = Path(path).expanduser()
        if not p.is_absolute() or not p.is_dir():
            raise ValueError('Select an existing folder using its full absolute path.')
        p = p.resolve()
        test = p / ('.songscope-write-check-' + uuid.uuid4().hex)
        with test.open('xb') as f:
            f.write(b'check')
        test.unlink()
        with self.lock:
            self.db['output'] = str(p)
            self.persist()
        return str(p)

    def capture_settings(self, data):
        camera = safe_text(data.get('camera', ''), 500)
        size = str(data.get('size', '1280x720'))
        if not re.fullmatch(r'\d{2,4}x\d{2,4}', size):
            raise ValueError('Use a resolution such as 1920x1080.')
        if not camera:
            raise ValueError('Select a camera.')
        fps = numeric(data.get('fps', 25), 1, 120, 'FPS')
        mode = data.get('input_format', 'mjpeg')
        if mode not in ('mjpeg', 'yuyv422'):
            raise ValueError('Unsupported camera input format.')
        camera_id = safe_text(data.get('camera_id', ''), 2000)
        saved = self.db.get('focus_profiles', {}).get(camera_id, {})
        focus = focus_request(data.get('focus', saved.get('focus', {'mode': 'device'})))
        return dict(camera=camera, camera_id=camera_id, size=size, fps=fps, input_format=mode, focus=focus,
                    transform=transform_settings(data.get('transform')))

    def focus_action(self, data):
        with self.command_lock:
            if self.active and self.active != 'preview':
                raise ValueError('Focus controls are locked during recording.')
            if self.synthetic:
                raise ValueError('Simulated capture cannot verify physical camera focus.')
            c = self.capture_settings(data)
            action = data.get('action')
            if action == 'save':
                if not self.project_file:
                    raise ValueError('Create or open a project before saving focus.')
                if self.active != 'preview' or self.live_capture != c or not (self.focus_status or {}).get('verified'):
                    raise ValueError('Apply and verify these settings in preview before saving.')
                with self.lock:
                    self.db.setdefault('focus_profiles', {})[c['camera_id']] = {
                        'camera': c['camera'], 'focus': copy.deepcopy(c['focus']), 'saved_at': now()}
                    self.persist()
                return copy.deepcopy(self.focus_status)
            if action not in ('read', 'apply'):
                raise ValueError('Unknown focus action.')
            self.stop_preview(); self.preempt_background()
            self.focus_status = None
            if action == 'read':
                bridge = self.focus_bridge_factory(ROOT, c)
                try:
                    self.focus_status = copy.deepcopy(bridge.report)
                    self.focus_status['verified'] = False
                    return copy.deepcopy(self.focus_status)
                finally:
                    bridge.close()
            if c['focus']['mode'] == 'device':
                raise ValueError('Choose Auto focus or Manual focus before applying.')
            self.launch(c, None, None)
            return copy.deepcopy(self.focus_status)

    def input_args(self, c):
        if self.synthetic:
            return ['-re', '-f', 'lavfi', '-i', f'testsrc2=size={c["size"]}:rate={c["fps"]}']
        fmt = ['-vcodec', 'mjpeg'] if c['input_format'] == 'mjpeg' else ['-pixel_format', 'yuyv422']
        return ['-f', 'dshow', '-rtbufsize', '256M', *fmt, '-video_size', c['size'],
                '-framerate', str(c['fps']), '-i', 'video=' + (c.get('camera_id') or c['camera'])]

    def stop_preview(self):
        if self.active == 'preview':
            self.request_stop('preview closed')
            if self.worker:
                self.worker.join(12)
            if self.proc:
                raise ValueError('The camera is still closing. Please try again shortly.')

    def preview_start(self, data):
        with self.command_lock:
            if self.active and self.active != 'preview':
                raise ValueError('Recording is in progress. Preview already uses the same video stream.')
            c = self.capture_settings(data)
            # Preview JPEGs stay uncropped. The browser adjusts framing immediately;
            # changing only framing must not reopen the camera or its focus bridge.
            with self.lock:
                if (self.active == 'preview' and self.proc and self.proc.poll() is None
                        and self.live_capture
                        and {k:v for k,v in self.live_capture.items() if k != 'transform'}
                            == {k:v for k,v in c.items() if k != 'transform'}):
                    self.live_capture = copy.deepcopy(c)
                    return
            self.stop_preview();self.preempt_background()
            self.launch(c, None, None)

    def start(self, data):
        with self.command_lock:
            if self.active and self.active != 'preview':
                raise ValueError('A recording or finalization is already in progress.')
            self.stop_preview();self.preempt_background()
            if not self.project_file:raise ValueError('Create or open a project before recording.')
            with self.lock:
                c = self.capture_settings(data)
                presets = json.loads((ROOT/'presets.json').read_text())
                assay = data.get('assay', 'OFT')
                if assay not in (*presets['seconds'], 'CUSTOM'):
                    raise ValueError('Choose a supported behavior.')
                custom_name = safe_text(data.get('custom_behavior',''),80).strip() if assay == 'CUSTOM' else ''
                if assay == 'CUSTOM' and not custom_name.strip(' .'):
                    raise ValueError('Enter a name for the custom behavioral test.')
                behavior_name = custom_name if assay == 'CUSTOM' else assay
                duration = numeric(data.get('duration_seconds', presets['seconds'].get(assay,360)),1,86400,'Duration (seconds)')
                test = data.get('test') is True
                if self.synthetic and not test: raise ValueError('Simulated capture requires test mode.')
                if not test and duration < 60: raise ValueError('Use test mode for durations under one minute.')
                output = Path(self.db['output'])
                if not self.db['output'] or not output.is_dir(): raise ValueError('Select an output folder first.')
                if shutil.disk_usage(output).free < duration * 4_000_000 + 1024**3:
                    raise ValueError('Insufficient space for MKV, MP4 and a 1 GB reserve.')
                prefix = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
                sequence = int(self.db.get('next_sequence',1))
                sid = f'{prefix}_{assay}_{sequence:03d}_{uuid.uuid4().hex[:8]}'
                self.db['next_sequence']=sequence+1
                mode = 'datetime_behavior'
                self.db['settings']=dict(capture=c,assay=assay,custom_behavior=custom_name,duration_seconds=duration,test=test,video_naming=mode)
                m = {k:safe_text(data.get(k,'')) for k in ('mouse_id','sex','cage','group','operator','notes')}
                if m['sex'] not in ('','F','M','Unknown'): raise ValueError('Invalid sex.')
                weight = data.get('weight_g')
                m['weight_g'] = numeric(weight,.001,1000,'Weight (g)') if weight not in (None,'') else None
                folder = output / ('TEST' if test else 'EXPERIMENT') / assay / sid
                folder.mkdir(parents=True,exist_ok=False)
                s = dict(id=sid, app_version=APP_VERSION, mouse=m, folder=str(folder),created_at=now(),test=test,
                         synthetic=self.synthetic, project_id=self.db['project']['id'], assay=assay,capture=c,recording_seconds=duration,
                         preset_source=presets, status='starting',phases={},events=[],phase='recording')
                auto_number = int(self.db.get('next_auto_video_id',1))
                automatic_id = f'Auto_ID{auto_number:02d}'
                relative = reserve_mp4(self.project_file.parent,self.db['sessions'],sid,prefix,m['mouse_id'],mode,test,behavior_name,automatic_id)
                if not m['mouse_id'].strip(' .'):
                    self.db['next_auto_video_id'] = auto_number+1
                s['behavior_name'] = behavior_name
                s['custom_behavior'] = custom_name
                s['video_naming'] = dict(mode=mode,mouse_id=m['mouse_id'],timestamp=prefix,behavior=behavior_name,automatic_id=automatic_id if not m['mouse_id'].strip(' .') else None)
                s['phases']['recording'] = dict(requested_at=now(),target_seconds=duration,file=sid+'.mkv',status='starting',mp4_relative=relative)
                self.db['sessions'].append(s)
                self.save_session(s); self.persist()
            self.launch(c,s,'recording')
            return sid

    def launch(self, c, session, phase):
        bridge = None
        self.focus_status = None
        self.live_capture = None
        try:
            if c.get('focus', {}).get('mode', 'device') != 'device':
                if self.synthetic:
                    raise ValueError('Simulated capture cannot verify physical camera focus.')
                bridge = self.focus_bridge_factory(ROOT, c, c['focus'])
                self.focus_status = copy.deepcopy(bridge.report)
                args = bridge.input_args()
            else:
                args = self.input_args(c)
            if session:
                session['focus'] = copy.deepcopy(self.focus_status) or {'verified': False, 'requested': c.get('focus'), 'reason': 'Camera default; focus not controlled.'}
                self.save_session(session)
        except Exception as exc:
            if bridge: bridge.close()
            self.focus_status = {'verified': False, 'camera_id': c.get('camera_id'), 'error': str(exc)}
            self.error = 'Recording blocked: ' + str(exc)
            if session:
                session.update(status='failed', error=self.error, focus=copy.deepcopy(self.focus_status))
                self.save_session(session); self.persist()
            raise ValueError(self.error) from exc
        cmd = [ffmpeg(), '-hide_banner', '-nostats', '-stats_period', '0.5', '-progress', 'pipe:2', *args]
        framing = framing_filter(c['size'], c.get('transform'))
        if session:
            duration = session[phase + '_seconds']
            # Preserve capture timestamps. No silent frame duplication to manufacture constant fps.
            cmd += ['-map', '0:v:0', '-an', '-t', str(duration), '-vf', framing, '-c:v', 'libx264', '-preset', 'veryfast',
                    '-crf', '18', '-maxrate', '16M', '-bufsize', '32M', '-pix_fmt', 'yuv420p',
                    '-fps_mode', 'passthrough', '-g', str(round(c['fps'] * 2)),
                    '-cluster_time_limit', '1000', '-flush_packets', '1', '-n', str(Path(session['folder']) / session['phases'][phase]['file'])]
        cmd += ['-map', '0:v:0', '-an']
        if session:
            cmd += ['-t', str(duration)]
        # Full source resolution while tuning focus; normal recording preview is
        # reduced so the existing encoder workload stays essentially unchanged.
        # Browser framing uses the same even-pixel crop geometry as the saved video.
        # Always send the full field, including during recording, to avoid double
        # cropping and allow instant adjustment without camera restarts.
        vf = 'fps=8' if not session else 'fps=8,scale=640:-2'
        cmd += ['-vf', vf, '-c:v', 'mjpeg', '-q:v', '6', '-threads', '1', '-f', 'image2pipe', 'pipe:1']
        with self.lock:
            self.active = session['id'] if session else 'preview'
            self.error = ''
            self.preview = b''
            self.preview_at = 0
            self.frames = 0
            self.media_seconds = 0
            self.started_mono = None
            self.stop_reason = None
            self.log_tail = []
            self.request_mono = self.last_progress = time.monotonic()
            try:
                self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE, creationflags=FLAGS, bufsize=0)
            except Exception as e:
                if bridge: bridge.close()
                self.active = None
                self.error = str(e)
                if session:
                    session['status'] = 'failed'
                    session['error'] = str(e)
                    self.save_session(session)
                    self.persist()
                raise
            self.live_capture = copy.deepcopy(c)
            self.worker = threading.Thread(target=self.monitor, args=(self.proc, session, phase, bridge), daemon=True)
            self.worker.start()

    def request_stop(self, reason):
        with self.lock:
            p = self.proc
            if p and p.poll() is None:
                self.stop_reason = reason
                try:
                    p.stdin.write(b'q\n')
                    p.stdin.flush()
                except (BrokenPipeError, OSError):
                    pass
                def kill_stuck():
                    try:
                        p.wait(8)
                    except subprocess.TimeoutExpired:
                        p.kill()
                threading.Thread(target=kill_stuck, daemon=True).start()

    def monitor(self, p, s, phase, bridge=None):
        def jpeg_reader():
            buf = b''
            while True:
                b = p.stdout.read(65536)  # bufsize=0: unbuffered pipe reads return available bytes
                if not b:
                    break
                buf += b
                while b'\xff\xd9' in buf:
                    end = buf.index(b'\xff\xd9') + 2
                    start = buf.find(b'\xff\xd8', 0, end)
                    if start >= 0:
                        self.preview = buf[start:end]
                        self.preview_at = time.monotonic()
                    buf = buf[end:]
                if len(buf) > 4_000_000:
                    buf = b''

        def stderr_reader():
            logfile = (Path(s['folder']) / (phase + '_capture.log')).open('wb') if s else None
            try:
                for raw in iter(p.stderr.readline, b''):
                    if logfile:
                        logfile.write(raw)
                        logfile.flush()
                    line = raw.decode('utf-8', 'replace').strip()
                    self.log_tail = (self.log_tail + [line])[-30:]
                    key, _, value = line.partition('=')
                    with self.lock:
                        if key == 'frame':
                            try:
                                self.frames = int(value)
                            except ValueError:
                                pass
                        if key == 'out_time_us':
                            try:
                                secs = max(0, int(value) / 1e6)
                            except ValueError:
                                continue
                            if secs > self.media_seconds:
                                self.last_progress = time.monotonic()
                            self.media_seconds = secs
                            if self.started_mono is None and secs > 0:
                                self.started_mono = time.monotonic() - secs
                                if s:
                                    # Host estimate; not a hardware exposure timestamp.
                                    stamp = dt.datetime.now().astimezone() - dt.timedelta(seconds=secs)
                                    s['phases'][phase]['first_frame_host_estimate'] = stamp.isoformat(timespec='milliseconds')
                                    s['phases'][phase]['status'] = s['status'] = 'recording'
                                    self.save_session(s)
                                    self.persist()
            finally:
                if logfile:
                    logfile.close()

        readers = [threading.Thread(target=jpeg_reader, daemon=True), threading.Thread(target=stderr_reader, daemon=True)]
        for t in readers:
            t.start()
        if os.name == 'nt':
            ctypes.windll.kernel32.SetThreadExecutionState(0x80000003)
        try:
            while p.poll() is None:
                time.sleep(0.5)
                if self.stop_reason:
                    continue
                elapsed = time.monotonic() - self.last_progress
                if elapsed > 20:
                    self.request_stop('No camera data or encoder progress for over 20 seconds')
                elif s and shutil.disk_usage(s['folder']).free < 256 * 1024**2:
                    self.request_stop('Free disk space is below 256 MB')
                elif s and self.started_mono and time.monotonic() - self.started_mono > s[phase + '_seconds'] + 15:
                    self.request_stop('Automatically stopped: target duration exceeded by more than 15 seconds')
            for t in readers:
                t.join(5)
            if bridge:
                bridge.close()
                if s:
                    s['focus_capture_log'] = bridge.logs[-30:]
            quality = dict(mp4_status='queued',qc='PENDING',qc_reasons=['Full timing QC will run with idle processing.'])
            with self.lock:
                if s:
                    record = s['phases'][phase]
                    path = Path(s['folder']) / record['file']
                    good = (p.returncode == 0 and self.stop_reason is None and self.frames > 0
                            and self.media_seconds >= s[phase + '_seconds'] - 0.6
                            and path.exists() and path.stat().st_size > 1024)
                    record.update(ended_at=now(), status='saved' if good else 'interrupted',
                                  frames=self.frames, media_seconds=self.media_seconds,
                                  bytes=path.stat().st_size if path.exists() else 0, returncode=p.returncode,
                                  average_fps=self.frames/self.media_seconds if self.media_seconds else None,
                                  reason=self.stop_reason,
                                  timestamp_note='Capture timestamps preserved; first-frame host time is estimated from progress, not hardware synchronization.')
                    expected = round(s[phase + '_seconds'] * s['capture']['fps'])
                    record['expected_frames'] = expected
                    record['frame_fraction'] = self.frames / expected if expected else 0
                    record['qc'] = 'REVIEW' if not good or abs(self.frames - expected) / max(expected,1) > .005 else 'frame_count_within_0.5_percent'
                    record['qc_note'] = 'Saving a file does not establish acquisition quality. A frame-count deviation over 0.5% requires review of exposure, lighting, USB and actual frame rate. Frame intervals also need review.'
                    record.update(quality)
                    if not good:
                        record['qc'] = 'REVIEW'
                    record['qc_note'] = 'Timing thresholds are screening criteria, not proof of analysis suitability. Missing-frame estimates are based on timestamps, not hardware counters. Check lighting, exposure, focus and USB.'
                    s['status'] = 'completed' if good else 'interrupted'
                    if not good:
                        s['error'] = self.stop_reason or '\n'.join(self.log_tail[-12:])
                        self.error = s['error']
                    self.save_session(s)
                    self.persist()
                elif p.returncode and not self.stop_reason:
                    self.error = '\n'.join(self.log_tail[-12:])
                self.proc = None
                self.active = None
                self.live_capture = None
                self.started_mono = None
                self.touch()
            if s and os.name == 'nt':
                import winsound
                winsound.MessageBeep(winsound.MB_OK if good else winsound.MB_ICONHAND)
        except Exception as e:
            self.error = 'Saving or monitoring failed: ' + str(e)
            self.request_stop(self.error)
            with self.lock:
                if s:
                    s['status'] = 'interrupted'
                    s['error'] = self.error
                    try:
                        self.persist()
                    except OSError:
                        pass
                self.active = None
                self.proc = None
        finally:
            if bridge and not bridge.closed.is_set():
                bridge.close()
            for pipe in (p.stdin, p.stdout, p.stderr):
                try:
                    pipe.close()
                except OSError:
                    pass
            if os.name == 'nt':
                ctypes.windll.kernel32.SetThreadExecutionState(0x80000000)

    def open_recording_video(self, data):
        with self.command_lock:
            if self.active and self.active != 'preview':
                raise ValueError('Wait for recording to finish before opening a video.')
            s = self.session(data['session_id'])
            q = s.get('phases', {}).get('recording', {})
            if q.get('mp4_status') != 'verified':
                raise ValueError('MP4 is not ready. Close preview and process pending videos first.')
            if q.get('mp4_relative'):
                if not self.project_file:
                    raise ValueError('Open the recording project first.')
                path = mp4_destination(self.project_file.parent, q['mp4_relative'])
            else:
                name = q.get('mp4_file')
                if not name or Path(name).name != name or Path(name).suffix.lower() != '.mp4':
                    raise ValueError('Invalid MP4 filename.')
                folder = Path(s['folder']).resolve()
                path = (folder/name).resolve()
                if path.parent != folder:
                    raise ValueError('MP4 path escapes the recording folder.')
            if not path.is_file():
                raise ValueError('MP4 file is missing. Restore the complete project folder.')
            os.startfile(str(path))
            return str(path)

    def manage_session(self, data):
        with self.command_lock:self.preempt_background()
        with self.lock:
            original = self.session(data['session_id'])
            if self.active == original['id'] or original['status'] not in ('completed', 'interrupted', 'failed'):
                raise ValueError('Finish or end the session before editing or deleting it.')
            s = copy.deepcopy(original)
            action = data.get('action')
            reason = safe_text(data.get('reason', ''), 1000)
            if not reason:
                raise ValueError('Enter a reason for this change.')
            before = dict(mouse=copy.deepcopy(s['mouse']),
                          deleted_at=s.get('deleted_at'))
            if action == 'edit':
                if s.get('deleted_at'):
                    raise ValueError('Restore the record before editing it.')
                fields = data.get('mouse', {})
                m = {k: safe_text(fields.get(k, s['mouse'].get(k, ''))) for k in
                     ('mouse_id', 'sex', 'cage', 'group', 'operator', 'notes')}
                if m['sex'] not in ('','F','M','Unknown'): raise ValueError('Invalid sex.')
                weight = fields.get('weight_g',s['mouse'].get('weight_g'))
                m['weight_g'] = numeric(weight,.001,1000,'Weight (g)') if weight not in (None,'') else None
                s['mouse'] = m
            elif action == 'delete':
                if s.get('deleted_at'): raise ValueError('This record is already deleted.')
                s['deleted_at'] = now()
            elif action == 'restore':
                if not s.get('deleted_at'): raise ValueError('This record is not deleted.')
                s.pop('deleted_at')
            else:
                raise ValueError('Unknown record action.')
            s.setdefault('revisions', []).append(dict(at=now(), action=action, reason=reason, before=before))
            s['updated_at'] = now()
            self.save_session(s)
            index = self.db['sessions'].index(original)
            self.db['sessions'][index] = s
            self.persist()
            return s['id']

    def event(self, data):
        with self.lock:
            s = self.session(data['session_id'])
            kind = data.get('kind', 'note')
            event = dict(at=now(), kind=kind, note=safe_text(data.get('note',''), 1000),
                         phase=s.get('phase'), video_seconds=self.media_seconds if self.active==s['id'] else None)
            if kind != 'note': raise ValueError('Only observation notes are supported.')
            s['events'].append(event)
            self.save_session(s)
            self.persist()


def export_log(db):
    out = io.StringIO(newline='')
    columns = ['record_id','assay','behavior_name','created_at','test','status','mouse_id','sex','cage','group','weight_g','operator','notes','duration_seconds','folder','qc','observed_fps','mp4_status','mp4_location']
    writer = csv.DictWriter(out,fieldnames=columns);writer.writeheader()
    for s in db['sessions']:
        if s.get('deleted_at'):continue
        q=s['phases'].get('recording',{})
        row=dict(record_id=s['id'],assay=s['assay'],behavior_name=s.get('behavior_name',s['assay']),created_at=s['created_at'],test=s['test'],status=s['status'],
                 **s['mouse'],duration_seconds=s['recording_seconds'],folder=s['folder'],
                 qc=q.get('qc',''),observed_fps=q.get('observed_fps',''),mp4_status=q.get('mp4_status',''),mp4_location=q.get('mp4_relative') or str(Path(s['folder'])/(q.get('mp4_file') or '')))
        # Prevent spreadsheet applications interpreting user text as formulas.
        row={k:("'"+v if isinstance(v,str) and v.startswith(('=','+','-','@')) else v) for k,v in row.items()}
        writer.writerow(row)
    return out.getvalue().encode('utf-8-sig')


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args): pass

    def reply(self, code, content, kind='application/json; charset=utf-8', attachment=None):
        if not isinstance(content,bytes):
            content=json.dumps(content,ensure_ascii=False,allow_nan=False).encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type',kind)
        self.send_header('Content-Length',str(len(content)))
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        if attachment: self.send_header('Content-Disposition',f'attachment; filename="{attachment}"')
        self.end_headers()
        try: self.wfile.write(content)
        except (BrokenPipeError,ConnectionResetError,ConnectionAbortedError): pass

    def allowed(self):
        port=self.server.server_port
        return self.headers.get('Host') in (f'127.0.0.1:{port}',f'localhost:{port}')

    def do_GET(self):
        if not self.allowed(): return self.reply(403,{'error':'Host rejected'})
        u=urlparse(self.path)
        if u.path=='/health':
            return self.reply(200,{'application':'SongScope','version':APP_VERSION})
        if u.path=='/favicon.ico':
            return self.reply(200,(ROOT/'assets/behavior-hub.ico').read_bytes(),'image/x-icon')
        if u.path=='/app-icon.png':
            return self.reply(200,(ROOT/'assets/behavior-hub.png').read_bytes(),'image/png')
        if u.path=='/':
            html=(ROOT/'index.html').read_text(encoding='utf-8').replace('__TOKEN__',TOKEN)
            html=html.replace('__FOCUS_UI__',(ROOT/'focus_ui.js').read_text(encoding='utf-8'))
            return self.reply(200,html.encode('utf-8'),'text/html; charset=utf-8')
        token=self.headers.get('X-SongScope-Token') or parse_qs(u.query).get('token',[''])[0]
        if not secrets.compare_digest(token,TOKEN): return self.reply(403,{'error':'Open the application from its local home page.'})
        try:
            if u.path=='/api/presets': return self.reply(200,json.loads((ROOT/'presets.json').read_text()))
            if u.path=='/api/state': return self.reply(200,self.server.rec.state())
            if u.path=='/api/devices':
                details = [{'name':'Synthetic test camera','id':'synthetic'}] if self.server.rec.synthetic else device_details()
                return self.reply(200,{'devices':[d['name'] for d in details], 'details':details})
            if u.path=='/api/preview':
                return self.reply(200,self.server.rec.preview,'image/jpeg') if self.server.rec.preview else self.reply(204,b'','image/jpeg')
            if u.path=='/api/export':
                return self.reply(200,export_log(self.server.rec.state()),'text/csv; charset=utf-8','BehaviorHub_record_log.csv')
            if u.path=='/api/json': return self.reply(200,self.server.rec.state(),attachment='BehaviorHub_records.json')
            return self.reply(404,{'error':'Not found'})
        except Exception as e: return self.reply(400,{'error':str(e)})

    def do_POST(self):
        if not self.allowed() or not secrets.compare_digest(self.headers.get('X-SongScope-Token',''),TOKEN):
            return self.reply(403,{'error':'Request rejected'})
        try:
            n=int(self.headers.get('Content-Length','0'))
            if n>100_000: raise ValueError('Request is too large.')
            data=json.loads(self.rfile.read(n) or b'{}')
            r=self.server.rec
            r.touch()
            route=urlparse(self.path).path
            result=None
            if route=='/api/project': result=r.project_action(data)
            elif route=='/api/queue': result=r.queue_action(data)
            elif route=='/api/activity': pass
            elif route=='/api/project-folder':
                import tkinter as tk
                from tkinter import filedialog
                root=tk.Tk();root.withdraw();root.attributes('-topmost',True)
                try: result=filedialog.askdirectory(title='Choose project parent folder',parent=root)
                finally: root.destroy()
            elif route=='/api/project-file':
                import tkinter as tk
                from tkinter import filedialog
                root=tk.Tk();root.withdraw();root.attributes('-topmost',True)
                try: result=filedialog.askopenfilename(title='Open Behavior Hub project',filetypes=[('Project files','*.project.json')],parent=root)
                finally: root.destroy()
            elif route=='/api/output': result=r.set_output(data['path'])
            elif route=='/api/folder':
                import tkinter as tk
                from tkinter import filedialog
                root=tk.Tk(); root.withdraw(); root.attributes('-topmost',True)
                try: path=filedialog.askdirectory(title='Select Behavior Hub recording folder',parent=root)
                finally: root.destroy()
                result=r.set_output(path) if path else None
            elif route=='/api/preview': r.preview_start(data)
            elif route=='/api/focus': result=r.focus_action(data)
            elif route=='/api/stop-preview':
                with r.command_lock: r.stop_preview()
            elif route=='/api/start': result=r.start(data)
            elif route=='/api/event': r.event(data)
            elif route=='/api/session': result=r.manage_session(data)
            elif route=='/api/open-video': result=r.open_recording_video(data)
            elif route=='/api/stop':
                reason=safe_text(data.get('reason',''))
                if not reason: raise ValueError('Enter a reason for stopping early.')
                r.request_stop(reason)
            elif route=='/api/options':
                if r.active: raise ValueError('Stop preview first. Camera modes cannot be queried during recording.')
                result=run_ff(['-list_options','true','-f','dshow','-i','video='+safe_text(data['camera'],500)]).stderr.decode('utf-8','replace')
            elif route=='/api/shutdown':
                if r.active and r.active!='preview': raise ValueError('Stop recording and wait for the file to be saved first.')
                r.stop_preview()
                r.close_services()
                threading.Thread(target=self.server.shutdown,daemon=True).start()
            else: return self.reply(404,{'error':'Not found'})
            return self.reply(200,{'ok':True,'result':result})
        except Exception as e:
            return self.reply(400,{'error':str(e) or type(e).__name__})

class LocalServer(ThreadingHTTPServer):
    # Windows SO_REUSEADDR can permit a second listener on an occupied port.
    # Exclusivity prevents the second launcher from briefly intercepting requests.
    allow_reuse_address = False
    def server_bind(self):
        if os.name == 'nt':
            self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        super().server_bind()


def existing_service(port):
    url = f'http://127.0.0.1:{port}'
    try:
        with urlopen(url+'/health',timeout=.6) as response:
            if json.load(response).get('application') == 'SongScope':
                return url
    except Exception:
        pass
    # Compatibility with the previous release, which did not expose /health.
    try:
        with urlopen(url,timeout=.6) as response:
            html=response.read(200_000).decode('utf-8','replace')
        if '<title>SongScope</title>' in html and 'X-SongScope-Token' in html:
            return url
    except Exception:
        pass
    return None


def show_startup_error(message, quiet=False):
    if quiet or os.name != 'nt':
        if sys.stderr:
            print(message,file=sys.stderr)
    else:
        ctypes.windll.user32.MessageBoxW(None,message,'Behavior Hub',0x10)
    return 1


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--port',type=int,default=43831)
    parser.add_argument('--data-dir',default=str(Path(os.environ.get('LOCALAPPDATA',str(Path.home())))/'SongScope'))
    parser.add_argument('--no-browser',action='store_true')
    parser.add_argument('--synthetic',action='store_true',help='TEST ONLY: generated frames, no camera')
    args=parser.parse_args()
    def handoff(port):
        url=existing_service(port)
        if url:
            if not args.no_browser:webbrowser.open(url)
            return True
        return False
    if handoff(args.port):
        return 0
    data_path=Path(args.data_dir)
    data_path.mkdir(parents=True,exist_ok=True)
    # Lock BEFORE binding or reading/writing the registry. A duplicate launch
    # must never recover/modify the live instance's sessions.
    data_lock=(data_path/'service.lock').open('a+b')
    instance_path=data_path/'instance.json'
    if os.name=='nt':
        import msvcrt
        try:
            data_lock.seek(0)
            msvcrt.locking(data_lock.fileno(),msvcrt.LK_NBLCK,1)
        except OSError:
            data_lock.close()
            # The first process may still be starting, or may use another port.
            for _ in range(30):
                ports=[args.port]
                try:
                    other=json.loads(instance_path.read_text(encoding='utf-8'))['port']
                    if isinstance(other,int) and 1<=other<=65535:ports.append(other)
                except (OSError,ValueError,KeyError):
                    pass
                if any(handoff(port) for port in dict.fromkeys(ports)):
                    return 0
                time.sleep(.1)
            return show_startup_error('The recording service is starting or shutting down. Please wait a few seconds and open Behavior Hub again. Your existing recordings have not been changed.',args.no_browser)
    server=None
    try:
        try:
            server=LocalServer(('127.0.0.1',args.port),Handler)
        except OSError:
            if handoff(args.port):return 0
            return show_startup_error(f'Port {args.port} is in use by another application. Close that application or start Behavior Hub with a different --port.',args.no_browser)
        server.rec=Recorder(args.data_dir,args.synthetic)
        atomic_json(instance_path,{'port':server.server_port,'pid':os.getpid(),'version':APP_VERSION})
        if not args.no_browser:webbrowser.open(f'http://127.0.0.1:{server.server_port}')
        server.serve_forever(poll_interval=.5)
    finally:
        if server:
            if hasattr(server,'rec'):
                server.rec.request_stop('Service shutdown')
                if server.rec.worker:server.rec.worker.join(12)
                server.rec.close_services()
            server.server_close()
        data_lock.close()
    return 0

if __name__=='__main__':
    try:
        code=main()
    except Exception as exc:
        code=show_startup_error('Behavior Hub could not start: '+str(exc),'--no-browser' in sys.argv)
    sys.exit(code)
