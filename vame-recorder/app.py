"""VAME Recorder: loopback-only Windows capture service. No cloud services."""
from background_jobs import BackgroundSupport
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

ROOT = Path(getattr(sys, '_MEIPASS', Path(__file__).parent))
FLAGS = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
TOKEN = secrets.token_urlsafe(32)
APP_VERSION = '1.3.0'

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
    r = run_ff(['-list_devices', 'true', '-f', 'dshow', '-i', 'dummy'])
    return re.findall(r'"([^"\r\n]+)" \(video\)', r.stderr.decode('utf-8', 'replace'))

def volume(weight):
    w = float(weight)
    if not math.isfinite(w) or w <= 0:
        raise ValueError('Weight must be a finite number greater than zero, in grams.')
    return w * 25 / 1000 / 2

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

class Recorder(BackgroundSupport):
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
        for s in self.db['sessions']:
            if s['status'] in ('starting', 'recording', 'finalizing'):
                s['status'] = 'interrupted'
                s['error'] = 'The service stopped before recording finished. Check the retained video files.'
                q=s.get('phases',{}).get(s.get('phase'))
                if q and s.get('app_version')=='1.3.0':q.update(status='interrupted',mp4_status='queued',qc='REVIEW')
                if Path(s['folder']).is_dir():self.save_session(s)
        self.ffmpeg_path=ffmpeg;self.process_flags=FLAGS
        for s in self.db['sessions']:
            for q in s['phases'].values():
                if q.get('mp4_status')=='processing':q['mp4_status']='queued'
        self.persist()
        self.init_background()

    def persist(self):
        atomic_json(self.db_path, self.db)

    def save_session(self, s):
        atomic_json(Path(s['folder']) / 'session.json', s)

    def session(self, sid):
        return next(s for s in self.db['sessions'] if s['id'] == sid)

    def state(self):
        with self.lock:
            d = copy.deepcopy(self.db)
            d.update(processing=self.queue_state(),preview_target_fps=self.preview_target_fps if self.active else 0)
            d.update(active=self.active, preview=bool(self.proc and self.active == 'preview'),
                     frames=self.frames, media_seconds=self.media_seconds,
                     elapsed=round(time.monotonic() - self.started_mono, 2) if self.started_mono else 0,
                     error=self.error, preview_age=round(time.monotonic() - self.preview_at, 1) if self.preview_at else None,
                     synthetic=self.synthetic)
            return d

    def mouse(self, data):
        with self.lock:
            m = {k: safe_text(data.get(k, '')) for k in ('mouse_id', 'sex', 'cage', 'group', 'operator', 'notes')}
            if not m['mouse_id'] or not m['cage'] or m['sex'] not in ('F', 'M', 'Unknown') or m['group'] not in ('KA', 'Saline'):
                raise ValueError('Enter a mouse ID, cage, sex and treatment group.')
            m['weight_g'] = numeric(data.get('weight_g'), 0.001, 1000, 'Weight (g)')
            m['volume_ml'] = volume(m['weight_g'])
            existing = next((x for x in self.db['mice'] if x['mouse_id'] == m['mouse_id']), None)
            m['updated_at'] = now()
            if existing:
                existing.update(m)
            else:
                if len(self.db['mice']) >= 12:
                    raise ValueError('This cohort supports up to 12 mice. Archive this cohort before using a new data directory.')
                self.db['mice'].append(m)
            self.persist()
            return m

    def set_output(self, path):
        p = Path(path).expanduser()
        if not p.is_absolute() or not p.is_dir():
            raise ValueError('Select an existing folder using its full absolute path.')
        p = p.resolve()
        test = p / ('.vame-write-check-' + uuid.uuid4().hex)
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
        return dict(camera=camera, size=size, fps=fps, input_format=mode)

    def input_args(self, c):
        if self.synthetic:
            return ['-re', '-f', 'lavfi', '-i', f'testsrc2=size={c["size"]}:rate={c["fps"]}']
        fmt = ['-vcodec', 'mjpeg'] if c['input_format'] == 'mjpeg' else ['-pixel_format', 'yuyv422']
        return ['-f', 'dshow', '-rtbufsize', '256M', *fmt, '-video_size', c['size'],
                '-framerate', str(c['fps']), '-i', 'video=' + c['camera']]

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
            self.stop_preview();self.preempt_background()
            self.launch(self.capture_settings(data), None, None)

    def start(self, data):
        with self.command_lock:
            if self.active and self.active != 'preview':
                raise ValueError('A recording is already in progress.')
            self.stop_preview();self.preempt_background()
            phase = data.get('phase')
            if phase not in ('baseline', 'post'):
                raise ValueError('Unknown recording phase.')
            with self.lock:
                if phase == 'baseline':
                    unfinished = [s for s in self.db['sessions'] if s['status'] in ('awaiting_injection', 'ready_post')]
                    if unfinished:
                        raise ValueError('Complete or end the pending session first.')
                    m = next(x for x in self.db['mice'] if x['mouse_id'] == data.get('mouse_id'))
                    c = self.capture_settings(data)
                    test = data.get('test') is True
                    if self.synthetic and not test:
                        raise ValueError('The simulated service can only create test sessions.')
                    baseline = numeric(data.get('baseline_seconds', 1500), 1, 86400, 'Baseline duration (seconds)')
                    post = numeric(data.get('post_seconds', 7200), 1, 86400, 'Post-return duration (seconds)')
                    if not test and (baseline < 60 or post < 60):
                        raise ValueError('Use test mode for durations under one minute.')
                    output = Path(self.db['output'])
                    if not self.db['output'] or not output.is_dir():
                        raise ValueError('Select an output folder first.')
                    # Budget for retained MKV AND MP4 at 16 Mbps each, plus 1 GiB reserve.
                    need = (baseline + post) * 4_000_000 + 1024**3
                    if shutil.disk_usage(output).free < need:
                        raise ValueError(f'Insufficient free space. At least {need/1024**3:.1f} GB is recommended.')
                    sid = dt.datetime.now().strftime('%Y%m%d_%H%M%S') + '_' + uuid.uuid4().hex[:8]
                    folder = output / ('TEST' if test else 'EXPERIMENT') / sid
                    folder.mkdir(parents=True, exist_ok=False)
                    s = dict(id=sid, app_version=APP_VERSION, mouse=copy.deepcopy(m), folder=str(folder), created_at=now(), test=test,
                             synthetic=self.synthetic, capture=c, baseline_seconds=baseline, post_seconds=post,
                             concentration_mg_ml=2, dose_mg_kg=25 if m['group']=='KA' else 0,
                             volume_factor_ml_g=0.0125, status='starting', phases={}, events=[], injection_at=None)
                    self.db['sessions'].append(s)
                else:
                    s = self.session(data['session_id'])
                    if s['status'] != 'ready_post':
                        raise ValueError('Confirm the injection before starting post-return recording.')
                    c = s['capture']
                    s['return_confirmed_at'] = now()
                duration = s[phase + '_seconds']
                s['phase'] = phase
                s['status'] = 'starting'
                s['phases'][phase] = {'requested_at': now(), 'target_seconds': duration,
                    'file': phase + '.mkv', 'status': 'starting'}
                self.save_session(s)
                self.persist()
            self.launch(c, s, phase)
            return s['id']

    def launch(self, c, session, phase):
        cmd = [ffmpeg(), '-hide_banner', '-nostats', '-stats_period', '0.5', '-progress', 'pipe:2', *self.input_args(c)]
        if session:
            duration = session[phase + '_seconds']
            # Preserve capture timestamps. No silent frame duplication to manufacture constant fps.
            cmd += ['-map', '0:v:0', '-an', '-t', str(duration), '-c:v', 'libx264', '-preset', 'veryfast',
                    '-crf', '18', '-maxrate', '16M', '-bufsize', '32M', '-pix_fmt', 'yuv420p',
                    '-fps_mode', 'passthrough', '-g', str(round(c['fps'] * 2)),
                    '-cluster_time_limit', '1000', '-flush_packets', '1', '-n', str(Path(session['folder']) / (phase + '.mkv'))]
        cmd += ['-map', '0:v:0', '-an']
        if session:
            cmd += ['-t', str(duration)]
        self.preview_target_fps=min(c['fps'],8 if session else 25)
        cmd += ['-vf', f'fps={self.preview_target_fps},scale=640:-2', '-c:v', 'mjpeg', '-q:v', '6', '-threads', '1', '-f', 'image2pipe', 'pipe:1']
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
                self.active = None
                self.error = str(e)
                if session:
                    session['status'] = 'failed'
                    session['error'] = str(e)
                    self.save_session(session)
                    self.persist()
                raise
            self.worker = threading.Thread(target=self.monitor, args=(self.proc, session, phase), daemon=True)
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

    def monitor(self, p, s, phase):
        def jpeg_reader():
            buf = b''
            while True:
                b = p.stdout.read(65536)
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
            quality = dict(mp4_status='queued',qc='PENDING',qc_reasons=['Full timing QC and MP4 conversion are queued for idle processing.'])
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
                    s['status'] = ('awaiting_injection' if phase == 'baseline' else 'completed') if good else 'interrupted'
                    if not good:
                        s['error'] = self.stop_reason or '\n'.join(self.log_tail[-12:])
                        self.error = s['error']
                    self.save_session(s)
                    self.persist()
                elif p.returncode and not self.stop_reason:
                    self.error = '\n'.join(self.log_tail[-12:])
                self.proc = None
                self.active = None
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
            for pipe in (p.stdin, p.stdout, p.stderr):
                try:
                    pipe.close()
                except OSError:
                    pass
            if os.name == 'nt':
                ctypes.windll.kernel32.SetThreadExecutionState(0x80000000)

    def manage_session(self, data):
        with self.command_lock:
            self.preempt_background()
            return self._manage_session(data)

    def _manage_session(self, data):
        with self.lock:
            original = self.session(data['session_id'])
            if self.active == original['id'] or original['status'] not in ('completed', 'interrupted', 'failed'):
                raise ValueError('Finish or end the session before editing or deleting it.')
            s = copy.deepcopy(original)
            action = data.get('action')
            reason = safe_text(data.get('reason', ''), 1000)
            if not reason:
                raise ValueError('Enter a reason for this change.')
            before = dict(mouse=copy.deepcopy(s['mouse']), actual_volume_ml=s.get('actual_volume_ml'),
                          deleted_at=s.get('deleted_at'))
            if action == 'edit':
                if s.get('deleted_at'):
                    raise ValueError('Restore the record before editing it.')
                fields = data.get('mouse', {})
                m = {k: safe_text(fields.get(k, s['mouse'].get(k, ''))) for k in
                     ('mouse_id', 'sex', 'cage', 'group', 'operator', 'notes')}
                if not any(x['mouse_id'] == m['mouse_id'] for x in self.db['mice']):
                    raise ValueError('Register the corrected mouse ID first, then select it here.')
                if not m['cage'] or m['sex'] not in ('F','M','Unknown') or m['group'] not in ('KA','Saline'):
                    raise ValueError('Enter a cage, valid sex and treatment group.')
                m['weight_g'] = numeric(fields.get('weight_g', s['mouse']['weight_g']), .001, 1000, 'Weight (g)')
                m['volume_ml'] = volume(m['weight_g'])
                actual = data.get('actual_volume_ml', s.get('actual_volume_ml'))
                if actual not in (None, ''):
                    s['actual_volume_ml'] = numeric(actual, .000001, 100, 'Actual volume (mL)')
                elif s.get('injection_at'):
                    raise ValueError('Actual volume is required for a logged injection.')
                else:
                    s.pop('actual_volume_ml', None)
                s['mouse'] = m
                s['dose_mg_kg'] = 25 if m['group'] == 'KA' else 0
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
            if kind == 'injection':
                if s['status'] != 'awaiting_injection':
                    raise ValueError('Injection cannot be confirmed again in the current state.')
                actual = numeric(data.get('actual_volume_ml'), 0.000001, 100, 'Actual injection volume (mL)')
                s['actual_volume_ml'] = actual
                s['injection_at'] = event['at']
                s['status'] = 'ready_post'
            elif kind == 'abandon':
                if self.active == s['id'] or s['status'] not in ('awaiting_injection','ready_post'):
                    raise ValueError('Only a waiting session can be ended here. Use Stop early during recording.')
                if not event['note']:
                    raise ValueError('Enter a reason for ending this session.')
                s['status'] = 'interrupted'
            elif kind != 'note':
                raise ValueError('Unknown event.')
            s['events'].append(event)
            self.save_session(s)
            self.persist()


def export_workbook(db):
    """Fill an artifact-authored template using OOXML, preserving layout/formulas."""
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    ET.register_namespace('', ns)
    def tag(s): return '{'+ns+'}'+s
    template = ROOT / 'VAME_Experimental_Log.xlsx'
    cells = {}
    caches = {}
    mice = db['mice'][:12]
    if len(db['mice']) > 12:
        raise ValueError('The print template supports 12 mice. Export JSON to retain all records.')
    for i,m in enumerate(mice):
        r = 8+i
        sessions = [s for s in db['sessions'] if s['mouse']['mouse_id']==m['mouse_id'] and not s['test'] and not s.get('deleted_at')]
        s = sessions[-1] if sessions else {}
        mm = s.get('mouse',m)
        vals = [mm['mouse_id'],mm['sex'],mm['cage'],mm['group'],mm['weight_g']]
        for col,val in zip('ABCDE',vals): cells[(1,f'{col}{r}')] = val
        caches[(1,f'F{r}')] = volume(mm['weight_g'])
        for col,val in [('G',s.get('actual_volume_ml','')),('H',s.get('created_at','')[:10]),('I',mm.get('operator','')),('J',s.get('status',''))]:
            cells[(1,f'{col}{r}')] = val
        sheet = i+2
        for cell,col in [('B5','A'),('F5','D'),('B6','B'),('F6','C'),('B7','E'),('F7','H'),('B8','F'),('F8','G'),('B9','I'),('F9','J')]:
            caches[(sheet,cell)] = caches.get((1,f'{col}{r}'),cells.get((1,f'{col}{r}'),''))
        for cell,val in [('B12',s.get('injection_at','')),('B13',s.get('return_confirmed_at','')),
                         ('B23',s.get('folder','')),('B24',s.get('id',''))]: cells[(sheet,cell)] = val
        if s:
            c=s['capture']
            cells[(sheet,'A3')]=f"Baseline: {s['baseline_seconds']/60:g} min; post-return: {s['post_seconds']/60:g} min; target {c['fps']:g} fps"
            cells[(sheet,'B21')]=f"{c['camera']}; {c['size']}; target {c['fps']} fps"
            cells[(sheet,'B14')]='; '.join(e.get('note','') for e in s.get('events',[]) if e['kind']=='injection')
        for phase,row in [('baseline',17),('post',19)]:
            p = s.get('phases',{}).get(phase,{})
            for cell,val in [(f'B{row}',p.get('first_frame_host_estimate','')),(f'F{row}',p.get('ended_at','')),
                             (f'B{row+1}',p.get('file','')),(f'F{row+1}',(p.get('status','')+' '+p.get('qc','')).strip())]: cells[(sheet,cell)]=val
        for j,e in enumerate(s.get('events',[])[:8]):
            cells[(sheet,f'A{28+j}')] = e['at'][11:19]
            cells[(sheet,f'B{28+j}')] = e.get('video_seconds')
            cells[(sheet,f'C{28+j}')] = e['kind']+': '+e.get('note','')
    out = io.BytesIO()
    with zipfile.ZipFile(template) as z, zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as dest:
        for name in z.namelist():
            raw = z.read(name)
            match = re.fullmatch(r'xl/worksheets/sheet(\d+)\.xml',name)
            if match:
                sn = int(match[1])
                changes = {cell:val for (sheet,cell),val in cells.items() if sheet==sn}
                formula_caches = {cell:val for (sheet,cell),val in caches.items() if sheet==sn}
                if changes or formula_caches:
                    root = ET.fromstring(raw)
                    data = root.find(tag('sheetData'))
                    for addr,val in changes.items():
                        rowno = re.search(r'\d+',addr).group()
                        row = next((x for x in data if x.get('r')==rowno),None)
                        if row is None: row = ET.SubElement(data,tag('row'),r=rowno)
                        c = next((x for x in row if x.get('r')==addr),None)
                        if c is None: c = ET.SubElement(row,tag('c'),r=addr)
                        for child in list(c): c.remove(child)
                        if isinstance(val,(float,int)):
                            c.attrib.pop('t',None)
                            ET.SubElement(c,tag('v')).text=str(val)
                        else:
                            c.set('t','inlineStr')
                            ET.SubElement(ET.SubElement(c,tag('is')),tag('t')).text='' if val is None else str(val)
                    for addr,val in formula_caches.items():
                        c=root.find('.//'+tag('c')+'[@r="'+addr+'"]')
                        if c is None: continue
                        for child in list(c):
                            if child.tag!=tag('f'):c.remove(child)
                        c.set('t','n' if isinstance(val,(float,int)) else 'str')
                        ET.SubElement(c,tag('v')).text=str(val)
                    raw=ET.tostring(root,encoding='utf-8',xml_declaration=True)
            dest.writestr(name,raw)
    return out.getvalue()


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
            return self.reply(200,{'application':'VAMERecorder','version':APP_VERSION})
        if u.path=='/':
            html=(ROOT/'index.html').read_text(encoding='utf-8').replace('__TOKEN__',TOKEN)
            return self.reply(200,html.encode('utf-8'),'text/html; charset=utf-8')
        token=self.headers.get('X-VAME-Token') or parse_qs(u.query).get('token',[''])[0]
        if not secrets.compare_digest(token,TOKEN): return self.reply(403,{'error':'Open the application from its local home page.'})
        try:
            if u.path=='/api/state': return self.reply(200,self.server.rec.state())
            if u.path=='/api/devices': return self.reply(200,{'devices':devices()})
            if u.path=='/api/preview':
                return self.reply(200,self.server.rec.preview,'image/jpeg') if self.server.rec.preview else self.reply(204,b'','image/jpeg')
            if u.path=='/api/export':
                return self.reply(200,export_workbook(self.server.rec.state()),'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet','VAME_Experimental_Log.xlsx')
            if u.path=='/api/json': return self.reply(200,self.server.rec.state(),attachment='VAME_registry.json')
            return self.reply(404,{'error':'Not found'})
        except Exception as e: return self.reply(400,{'error':str(e)})

    def do_POST(self):
        if not self.allowed() or not secrets.compare_digest(self.headers.get('X-VAME-Token',''),TOKEN):
            return self.reply(403,{'error':'Request rejected'})
        try:
            n=int(self.headers.get('Content-Length','0'))
            if n>100_000: raise ValueError('Request is too large.')
            data=json.loads(self.rfile.read(n) or b'{}')
            r=self.server.rec
            r.touch()
            route=urlparse(self.path).path
            result=None
            if route=='/api/queue': result=r.queue_action(data)
            elif route=='/api/activity': pass
            elif route=='/api/mouse': result=r.mouse(data)
            elif route=='/api/output': result=r.set_output(data['path'])
            elif route=='/api/folder':
                import tkinter as tk
                from tkinter import filedialog
                root=tk.Tk(); root.withdraw(); root.attributes('-topmost',True)
                try: path=filedialog.askdirectory(title='Select VAME recording folder',parent=root)
                finally: root.destroy()
                result=r.set_output(path) if path else None
            elif route=='/api/preview': r.preview_start(data)
            elif route=='/api/stop-preview':
                with r.command_lock: r.stop_preview()
            elif route=='/api/start': result=r.start(data)
            elif route=='/api/event': r.event(data)
            elif route=='/api/session': result=r.manage_session(data)
            elif route=='/api/stop':
                reason=safe_text(data.get('reason',''))
                if not reason: raise ValueError('Enter a reason for stopping early.')
                r.request_stop(reason)
            elif route=='/api/options':
                if r.active: raise ValueError('Stop preview first. Camera modes cannot be queried during recording.')
                result=run_ff(['-list_options','true','-f','dshow','-i','video='+safe_text(data['camera'],500)]).stderr.decode('utf-8','replace')
            elif route=='/api/shutdown':
                if r.active and r.active!='preview': raise ValueError('Stop recording and wait for the file to be saved first.')
                r.stop_preview();r.close_services()
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
            if json.load(response).get('application') == 'VAMERecorder':
                return url
    except Exception:
        pass
    # Compatibility with the previous release, which did not expose /health.
    try:
        with urlopen(url,timeout=.6) as response:
            html=response.read(200_000).decode('utf-8','replace')
        if '<title>VAME Recorder</title>' in html and 'X-VAME-Token' in html:
            return url
    except Exception:
        pass
    return None


def show_startup_error(message, quiet=False):
    if quiet or os.name != 'nt':
        if sys.stderr:
            print(message,file=sys.stderr)
    else:
        ctypes.windll.user32.MessageBoxW(None,message,'VAME Recorder',0x10)
    return 1


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--port',type=int,default=43821)
    parser.add_argument('--data-dir',default=str(Path(os.environ.get('LOCALAPPDATA',str(Path.home())))/'VAMERecorder'))
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
            return show_startup_error('The recording service is starting or shutting down. Please wait a few seconds and open VAMERecorder again. Your existing recordings have not been changed.',args.no_browser)
    server=None
    try:
        try:
            server=LocalServer(('127.0.0.1',args.port),Handler)
        except OSError:
            if handoff(args.port):return 0
            return show_startup_error(f'Port {args.port} is in use by another application. Close that application or start VAME Recorder with a different --port.',args.no_browser)
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
        code=show_startup_error('VAME Recorder could not start: '+str(exc),'--no-browser' in sys.argv)
    sys.exit(code)
