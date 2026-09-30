"""Reproducible minimal DirectShow -> NUT bridge. Run on Windows.

Downloads pinned upstream source/toolchain into ignored build/, never installs
system software. Only the capture bridge changes; the encoder stays unchanged.
"""
from pathlib import Path
import hashlib
import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / 'build/focus-tools'
VERSION = '7.1.3'

def build():
    WORK.mkdir(parents=True, exist_ok=True)
    archive = WORK / f'ffmpeg-{VERSION}.tar.xz'
    if not archive.exists():
        urllib.request.urlretrieve(f'https://ffmpeg.org/releases/ffmpeg-{VERSION}.tar.xz', archive)
    tool = WORK / 'w64devkit.7z.exe'
    if not tool.exists():
        urllib.request.urlretrieve('https://github.com/skeeto/w64devkit/releases/download/v2.10.0/w64devkit-x64-2.10.0.7z.exe', tool)
    for path, expected in ((archive, 'f0bf043299db9e3caacb435a712fc541fbb07df613c4b893e8b77e67baf3adbe'),
                           (tool, '18d0a4c71a166f8401ab6305781bec5882b40b5e06ba9807c61cb5f3b3c6325e')):
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError('Build input hash mismatch: ' + str(path))
    tools = WORK / 'w64devkit/bin'
    if not (tools / 'gcc.exe').exists():
        subprocess.run([str(tool), '-y', '-o' + str(WORK)], check=True)
    # FFmpeg's makefiles require a source/build path without spaces.
    build_root = Path(tempfile.gettempdir()) / 'BehaviorHub-focus-build'
    build_root.mkdir(exist_ok=True)
    source = build_root / f'ffmpeg-{VERSION}'
    with tarfile.open(archive) as tar:
        # Restore only patched files on rebuild; retain compiled objects.
        for name in ('libavdevice/dshow.c', 'libavdevice/dshow_capture.h'):
            member = tar.getmember(f'ffmpeg-{VERSION}/{name}')
            if not source.exists():
                tar.extractall(build_root, filter='data')
                break
            (source / name).write_bytes(tar.extractfile(member).read())
    header = source / 'libavdevice/dshow_capture.h'
    text = header.read_text()
    text = text.replace('    int   use_video_device_timestamps;',
                        '    int   use_video_device_timestamps;\n    int bh_focus_mode;\n    int bh_focus_value;\n    volatile LONG bh_focus_ready;')
    header.write_text(text)
    c = source / 'libavdevice/dshow.c'
    text = c.read_text()
    insert = (ROOT / 'native/focus_control.c').read_text()
    text = text.replace('static int\ndshow_read_close', insert + '\nstatic int\ndshow_read_close', 1)
    assert insert in text
    text = text.replace('//    dump_videohdr(s, vdhdr);',
                        'if (ctx->bh_focus_mode >= 0 && !InterlockedCompareExchange(&ctx->bh_focus_ready, 0, 0))\n        return;')
    # Apply and verify AFTER Run, which may reset device settings.
    old = '    ret = 0;\n\nerror:\n\n    if (devenum)'
    assert old in text
    text = text.replace(old, '    ret = bh_focus_setup(avctx);\n\nerror:\n\n    if (devenum)', 1)
    text = text.replace('static const AVOption options[] = {', '''static const AVOption options[] = {
    { "bh_focus_mode", "Behavior Hub focus: -1 unmanaged, 0 read, 1 auto, 2 manual", OFFSET(bh_focus_mode), AV_OPT_TYPE_INT, {.i64 = -1}, -1, 2, DEC },
    { "bh_focus_value", "Device focus value", OFFSET(bh_focus_value), AV_OPT_TYPE_INT, {.i64 = 0}, INT_MIN, INT_MAX, DEC },''')
    c.write_text(text)
    env = os.environ.copy()
    env['PATH'] = str(tools) + os.pathsep + env['PATH']
    env['TMPDIR'] = build_root.as_posix()
    configure = ['./configure', '--target-os=mingw32', '--arch=x86_64', '--cc=gcc', '--disable-everything', '--disable-autodetect', '--disable-x86asm',
                 '--disable-doc', '--disable-debug', '--disable-network', '--disable-ffplay', '--disable-ffprobe',
                 '--enable-indev=dshow', '--enable-muxer=nut', '--enable-protocol=pipe',
                 '--enable-decoder=mjpeg,rawvideo', '--enable-parser=mjpeg',
                 '--enable-filter=null', '--extra-version=BehaviorHub-focus1']
    if not (source / 'ffbuild/config.mak').exists():
        subprocess.run([str(tools / 'sh.exe'), *configure], cwd=source, env=env, check=True)
    subprocess.run([str(tools / 'make.exe'), '-j8', 'V=1'], cwd=source, env=env, check=True)
    shutil.copy2(source / 'ffmpeg.exe', ROOT / 'focus_capture.exe')
    shutil.copy2(archive, ROOT / 'native' / archive.name)
    shutil.copy2(WORK / 'w64devkit/COPYING.MinGW-w64-runtime.txt', ROOT / 'native/mingw-w64-COPYING')
    for name in ('COPYING.LGPLv2.1', 'LICENSE.md'):
        shutil.copy2(source / name, ROOT / 'native' / name)
    print('FOCUS_BRIDGE_BUILT', hashlib.sha256((ROOT / 'focus_capture.exe').read_bytes()).hexdigest())

if __name__ == '__main__':
    build()
