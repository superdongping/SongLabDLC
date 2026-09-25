# VAME Recorder 1.3.0

A local Windows application for bottom-view mouse video acquisition and printable experiment logs. No account or cloud service is required.

## Install or upgrade

1. Extract `VAMERecorder_Windows_English_1.3.0.zip` into a new local folder. Keep `VAMERecorder.exe` and the entire `_internal` folder together.
2. When upgrading, finish any recording in the previous version, then click **Exit application** in that version. Closing the browser tab alone leaves its recording service running.
3. Double-click the new `VAMERecorder.exe`. It opens `http://127.0.0.1:43821` in your default browser. Python and FFmpeg do not need to be installed separately.
4. If the service is already running, double-clicking the executable simply reopens it. The existing session and recordings are not changed.
5. On a different Windows computer, copy the entire extracted folder, connect the camera, and select that computer's output folder. Supported platform: Windows 10/11 x64.

Existing mouse records and settings remain in `%LOCALAPPDATA%\VAMERecorder`. The upgrade does not translate, delete or modify your previously entered notes, mouse IDs or saved videos. If the old Chinese page still appears, exit the old service first and then launch the new executable.

## Recording workflow

1. Enter **Mouse ID**, **Cage ID**, **Sex**, **Group**, **Weight (g)** and optionally the operator and notes. Click **Save mouse**. A cohort supports up to 12 mice, with a planned allocation of 6 KA and 6 saline. The program does not assign treatment groups automatically.
2. Verify the calculated injection volume. KA is fixed at **2 mg/mL** and **25 mg/kg**. Volume in mL = body weight in grams × **0.0125**. Saline uses the same volume-per-weight factor. For example, 25 g gives 0.3125 mL (312.5 microliters).
3. Choose an existing output folder. Select the camera, resolution, frame rate and input format. Defaults are 1280×720, 25 fps, MJPEG. Audio is not recorded.
4. Click **Open preview** and check the field of view, lighting and focus. Click **Start baseline**. The default baseline lasts 25 minutes and stops automatically, saving `baseline.mkv`.
5. The application displays the appropriate KA or saline injection prompt. Enter the **actual injection volume**, add any lot/route notes, and click **Injection complete: record time** immediately after injection. This records the current computer time; the program does not administer treatment.
6. Return the mouse to the arena and click **Mouse returned: start second phase**. A full 120 minutes is recorded by default, beginning with this phase rather than the earlier injection time. It stops automatically and saves `post.mkv`.
7. Use **Log event time** for observations, interventions or deviations. To stop early, click **Stop early and save** and enter a reason. The session is marked interrupted and existing video is retained.
8. Complete the current session, or enter a reason and click **End waiting session**, before beginning the next mouse.

Refreshing or closing the web page does not stop a running recording. Reopen the executable to reconnect. To stop the background service, use **Exit application** after recording has stopped and saving has finished.

## Test mode

Select **Test mode** to use seconds instead of minutes. Defaults change to a 10-second baseline and a 15-second second phase. Both values are editable. Test sessions are saved under `TEST`; non-test sessions are saved under `EXPERIMENT`. Test injection entries are simulated software events and should be labeled accordingly in the notes.

## Files and storage

Each session receives a unique timestamp-and-random-ID folder. Files are not overwritten and treatment groups are not included in video filenames.

- `baseline.mkv` and `post.mkv`: H.264 video in an MKV container, written continuously to disk.
- `session.json`: mouse information snapshot, calculation parameters, actual volume, injection/return times, recording information, events and quality flags.
- `baseline_capture.log` and `post_capture.log`: camera negotiation, encoder progress and diagnostic messages.

Mouse records, the session index and the selected output folder are stored in `%LOCALAPPDATA%\VAMERecorder\registry.json`. Back up this file and the actual video folders. The **Download all records (JSON)** button exports the full index. Moving files to another computer does not automatically update their stored paths.

The service runs only on `127.0.0.1`; other computers cannot connect remotely. Each acquisition computer runs its own copy.

## Printable Excel logs

**Download printable Excel** generates a cohort table and one sheet per mouse, up to 12 mice. Unknown information stays blank. Each mouse page uses the latest non-test session; if none exists, it uses the current registered mouse information. Test sessions do not populate experiment logs.

The workbook uses Letter paper, with a landscape cohort table and one portrait page per mouse. In Excel, select **Print Entire Workbook**. Yellow cells in the cohort table accept input; calculated volumes and mouse pages update through formulas. Do not sort cohort rows independently because the mouse pages refer to fixed rows.

Each printed page displays up to eight events. The full event history remains in JSON. The standalone blank workbook is included for printing before entering mouse information.

## Timing and acquisition quality

**Saved does not mean acquisition quality has passed.** The software compares encoded frame counts against target duration × target fps. A deviation greater than 0.5% is flagged **REVIEW**. It preserves capture timestamps and does not silently duplicate frames to manufacture the target rate. Frame intervals, image clarity, occlusion and decoding must also be checked before analysis.

The displayed countdown follows encoder progress and may lag capture. First-frame computer time is estimated from progress and buffering; it is not a hardware-synchronized exposure timestamp. Injection and return button times are recorded separately with the local UTC offset.

The original C920 short tests produced fewer frames than the requested 25 fps, including timestamp gaps up to 120 ms. Verify the actual camera, USB connection, lighting and exposure settings on the acquisition computer before claiming constant 25 fps. See `VALIDATION.md` for the scope and limits of testing.

## Interruptions

- No camera/encoder progress for about 20 seconds, low disk space, an encoder failure or manual early stop marks the session interrupted. Existing files are retained.
- On restart, an unfinished session is marked interrupted. Check the retained video; abrupt interruption may require container-index repair and full recovery is not guaranteed.
- Pause, seamless resume and automatic repeat recording are not supported. Use a new session for a repeat; retain the original session.
- The service requests that Windows not automatically sleep during capture. This cannot prevent manual sleep, shutdown, unplugging or power loss.

## Scope and parameters

Study parameters were confirmed by the project owner on 2026-09-21: 6 KA and 6 saline mice; bottom view only; no EEG or side camera; approximately 25-minute baseline, then injection outside the arena, followed by 120-minute post-return video; KA 2 mg/mL at 25 mg/kg; volume-matched saline. The program records the experiment and does not perform pose estimation, VAME analysis, Racine scoring or electrographic seizure detection.

## Developer information

Application source is included under `Source`. Development uses Python 3.11, imageio-ffmpeg 0.6.0 and PyInstaller 6.22.3. Copy `_internal/ffmpeg.exe` into `Source` to use that binary when running the Python source, or install imageio-ffmpeg. `test_recorder.py` exercises the state machine; `test_http.py` tests the packaged executable against an isolated simulated camera.

Third-party component information is in `THIRD_PARTY_NOTICES.txt` and `licenses/`.

## Version 1.3.0: MP4 and timing verification

Each phase retains its original H.264 MKV and queues MP4 stream copy (no re-encoding) plus full timing QC. The injection prompt, second phase and next animal do not wait for conversion. This can take several minutes for long recordings. Both files are retained, requiring approximately twice the video storage. Old recordings are not automatically converted.

The MP4 is named baseline.mp4 or post.mp4 only after full decoding confirms identical frame hashes/counts and relative timestamps within 1.1 ms of MKV. Failed conversions retain the MKV and may leave a .pending.mp4 that must not be treated as verified. An interrupted acquisition remains interrupted even if its retained video converts successfully.

Per-phase *_frame_timing.csv includes frame index (zero-based), relative presentation timestamp, interval and decoded-frame MD5. Session JSON records measured FPS as (frames - 1) / timestamp span, maximum interval, intervals over 1.5 times the target period, estimated missing frames from gaps and target frame-count deficit. These estimates are not hardware drop counters. Count deviation above 0.5%, long gaps, non-increasing timestamps or processing errors trigger REVIEW. TIMING_CHECKS_PASSED is only a timing screen, not proof of image clarity or suitability for analysis. At 25 fps the nominal interval is 40 ms and the long-gap threshold is 60 ms.

Validate lighting, exposure, focus and USB bandwidth on the actual recording computer. Changing the container does not fix dropped frames; the application never fabricates frames to meet the requested frame rate.

## Record management (1.3.0)

New recordings default to 1280 x 720. Existing session capture settings are unchanged.

Use Edit on a completed, interrupted or failed session to correct mouse ID (choose an already registered ID), sex, cage, group, weight, operator, notes and actual injection volume. Register a new mouse ID first if needed. Edits affect this session only, not the roster or other sessions. Calculated volume updates from weight; actual administered volume is not silently recalculated. Video files, capture timestamps and quality results are immutable through this editor. A correction reason is required and prior values are retained in session JSON revisions.

Delete removes a record from the normal list and excludes it from printable Excel. It does NOT delete the mouse registration, MKV, MP4, timing CSV or session folder. Select Show deleted records and Restore to undo. JSON retains all records including deletion metadata. Finish or end a session before editing/deleting it; active and waiting sessions are protected.

## Version 1.3.0: idle queue and smoother preview

- Setup-only preview targets up to 25 fps (limited by the selected camera FPS), replacing the old 3 fps preview. During recording, preview targets up to 8 fps at width 640 to reduce load. These are preview targets, not measured camera acquisition FPS. Video recording settings and timestamp preservation are unchanged.
- Automatic MP4 and full frame-QC work starts only with the camera closed, no pending injection/second-phase session, and no recent application interaction for the idle delay (default 60 seconds; adjustable 5-3600 seconds). This is application idle detection, not a system-wide CPU-idle measurement.
- Idle processing queue offers Process pending videos, Pause processing, Resume automatic and Retry failed. Manual processing can be requested between phases after closing preview. Starting recording or preview cancels any active processing job, returns it to the queue, and restarts it later. This is job restart, not byte-level resume.
- Capture is complete when its MKV closes. MP4: queued and QC: PENDING do not mean acquisition quality has passed. Check the final timing report after processing. MKV is never deleted by conversion.
- The queue persists across service restarts. Exit application safely stops processing and leaves pending work for reopening. Closing the browser alone leaves the service running. Deleted records are excluded from automatic processing. Old recordings without queued jobs are not automatically converted.
- Keep the service open after the experiment with preview closed, or manually process pending videos. Conversion failures do not block the next animal. Resolve disk space/missing-file issues and use Retry failed.

Validation uses synthetic videos and isolated state; actual camera FPS, image clarity and full-duration performance still require testing on the recording computer. Preview smoothness does not establish recording quality.

## Source and Windows build

This repository contains source and an unfilled printable log template. Experimental recordings, mouse registries and local outputs are not included. VAME Recorder uses port 43821 and its own local registry; Behavior Hub remains separate on port 43831.

From this directory, with Python 3.11 on Windows:

```powershell
python -m pip install -r requirements.txt
python app.py
```

To build the portable executable and ZIP:

```powershell
python build_windows.py
```

This creates `dist/1.3.0/VAMERecorder/VAMERecorder.exe` and `VAMERecorder_Windows_English_1.3.0.zip`. FFmpeg is supplied by imageio-ffmpeg. Keep the complete extracted folder together.

Synthetic regression tests (no physical camera):

```powershell
python -m unittest test_recorder test_media_quality test_idle -v
npm ci --ignore-scripts
python test_http.py
python test_package_dom.py
python test_relaunch.py
```

The packaged tests require a completed build and use isolated temporary state. DOM emulation is not real-browser visual validation. Review VALIDATION.md for remaining hardware and long-duration checks.
