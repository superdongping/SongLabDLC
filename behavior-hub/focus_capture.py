"""Focus-gated native capture, lossless NUT transport to the existing encoder.

The native bridge sets/reads IAMCameraControl on its active DirectShow source
filter AFTER starting the graph. It withholds all video until readback succeeds.
The bridge stream-copies camera packets and their timestamps; no synthetic FPS.
"""
import copy
import datetime as dt
import json
import os
from pathlib import Path
import socket
import subprocess
import threading

FLAGS = subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0


def focus_request(value):
    if not isinstance(value, dict):
        raise ValueError('Invalid focus settings.')
    mode = value.get('mode')
    if mode not in ('auto', 'manual', 'device'):
        raise ValueError('Choose automatic or manual focus.')
    result = {'mode': mode}
    if mode == 'manual':
        v = value.get('value')
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not -2147483648 <= v <= 2147483647 or int(v) != v:
            raise ValueError('Manual focus must be an integer reported by the device.')
        result['value'] = int(v)
    return result


def validate_report(report, request):
    if not isinstance(report, dict) or report.get('protocol') != 1 or report.get('source') != 'active_capture_filter':
        raise ValueError('Cannot verify the active camera connection. Recording blocked.')
    for key in ('minimum', 'maximum', 'step', 'value', 'flags', 'capabilities'):
        if type(report.get(key)) is not int:
            raise ValueError('Incomplete focus readback. Recording blocked.')
    if report['step'] <= 0 or report['maximum'] < report['minimum']:
        raise ValueError('Invalid device focus range.')
    if request is not None:
        flag = 1 if request['mode'] == 'auto' else 2
        if report.get('verified') is not True or report['flags'] & 3 != flag or not report['capabilities'] & flag:
            raise ValueError('Focus mode readback failed. Recording blocked; apply focus again.')
        if flag == 2:
            value = request['value']
            if not report['minimum'] <= value <= report['maximum'] or (value-report['minimum']) % report['step'] or value != report['value']:
                raise ValueError('Manual focus readback failed. Recording blocked; apply focus again.')
    return report


class FocusBridge:
    def __init__(self, root, capture, request=None):
        self.request = copy.deepcopy(request)
        self.capture = copy.deepcopy(capture)
        self.report = None
        self.error = ''
        self.logs = []
        self.ready = threading.Event()
        self.closed = threading.Event()
        self.listener = None
        self.connection = None
        self.relay_thread = None
        exe = Path(root) / 'focus_capture.exe'
        if not exe.is_file():
            raise ValueError('Focus capture component is missing. Extract the complete application package.')
        camera_id = capture.get('camera_id', '')
        if not camera_id.startswith('@device_'):
            raise ValueError('Refresh cameras and select the physical camera before controlling focus.')
        mode = 0 if request is None else 1 if request['mode'] == 'auto' else 2
        fmt = ['-vcodec', 'mjpeg'] if capture['input_format'] == 'mjpeg' else ['-pixel_format', 'yuyv422']
        cmd = [str(exe), '-hide_banner', '-nostdin', '-nostats', '-f', 'dshow', '-rtbufsize', '256M',
               *fmt, '-video_size', capture['size'], '-framerate', str(capture['fps']),
               '-bh_focus_mode', str(mode), '-bh_focus_value', str((request or {}).get('value', 0)),
               '-i', 'video=' + camera_id, '-map', '0:v:0', '-an', '-c:v', 'copy',
               '-f', 'nut', '-flush_packets', '1', 'pipe:1']
        self.proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, creationflags=FLAGS, bufsize=0)
        self.reader = threading.Thread(target=self._read_status, daemon=True)
        self.reader.start()
        try:
            if not self.ready.wait(15):
                raise ValueError('Camera focus verification timed out. Recording blocked.')
            if self.report is None or self.proc.poll() is not None:
                raise ValueError(self.error or 'Camera focus control unavailable. Close other camera applications and try again. ' + '\n'.join(self.logs[-4:]))
            validate_report(self.report, request)
            self.report.update(camera_id=camera_id, requested=self.request,
                               checked_at=dt.datetime.now().astimezone().isoformat(timespec='milliseconds'))
        except Exception:
            self.close()
            raise

    def _read_status(self):
        try:
            for raw in iter(self.proc.stderr.readline, b''):
                line = raw.decode('utf-8', 'replace').strip()
                self.logs.append(line)
                self.logs = self.logs[-60:]
                if 'BH_FOCUS_ERROR ' in line:
                    self.error = line.split('BH_FOCUS_ERROR ', 1)[1]
                    self.ready.set()
                elif 'BH_FOCUS {' in line:
                    try:
                        self.report = json.loads(line.split('BH_FOCUS ', 1)[1])
                    except ValueError:
                        self.error = 'Invalid camera focus response.'
                    self.ready.set()
        finally:
            self.ready.set()

    def input_args(self):
        # One loopback-only connection transports unmodified NUT packets. Encoder
        # stdin remains free for FFmpeg's graceful q/stop command.
        self.listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if os.name == 'nt':
            self.listener.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        self.listener.bind(('127.0.0.1', 0))
        self.listener.listen(1)
        self.listener.settimeout(15)
        port = self.listener.getsockname()[1]
        self.relay_thread = threading.Thread(target=self._relay, daemon=True)
        self.relay_thread.start()
        return ['-f', 'nut', '-i', f'tcp://127.0.0.1:{port}']

    def _relay(self):
        try:
            self.connection, _ = self.listener.accept()
            self.listener.close()
            self.connection.settimeout(10)
            while not self.closed.is_set():
                data = self.proc.stdout.read(65536)
                if not data:
                    break
                self.connection.sendall(data)
        except (OSError, ValueError) as exc:
            if not self.closed.is_set():
                self.error = 'Camera transport closed: ' + str(exc)
        finally:
            if self.connection:
                self.connection.close()

    def close(self):
        self.closed.set()
        for sock in (self.listener, self.connection):
            if sock:
                try:
                    sock.shutdown(socket.SHUT_RDWR)
                except OSError:
                    pass
                sock.close()
        if self.proc.poll() is None:
            self.proc.terminate()
        self.proc.wait(5)
        self.reader.join(3)
        if self.relay_thread:
            self.relay_thread.join(3)
        for pipe in (self.proc.stdout, self.proc.stderr):
            pipe.close()
