# Behavior Hub 1.1.0 — Behavioral Recording Studio

Behavior Hub is an English-language, local Windows application for behavioral video capture. It is independent of VAME Recorder and does not administer treatment or launch DLC/MATLAB analysis.

## Start

1. Extract the complete ZIP into a local folder. Keep BehaviorHub.exe and _internal together.
2. Open BehaviorHub.exe. The browser opens http://127.0.0.1:43831. Reopening the EXE reconnects to the running service.
3. Create a project (name and parent folder), or open BehaviorHub.project.json from an existing project. Choose a camera and behavior. Default capture is 1280x720, target 25 fps, MJPEG, no audio.
4. Confirm recording duration, check framing/light/focus with Open preview, then Start recording. All mouse fields are optional.
5. Recording stops automatically. The MKV closes and the record is queued for background processing. You can start the next animal without waiting for MP4 or full frame QC.
6. Exit application stops the background service. Closing a browser tab alone does not stop capture.

State is in %LOCALAPPDATA%\SongScope. VAME Recorder state and existing experiments are not modified. Do not try to open the same physical webcam in both applications at once.

## Presets

| Behavior | Default recording preset |
|---|---:|
| OFT | 6 minutes |
| NPR | 6 minutes |
| Zero Maze | 6 minutes |
| Y-maze | 8 minutes |
| FST | 5 minutes |
| TST | 6 minutes |

Presets were imported from SongLabDLC helpers/get_default_behavior_options.m at commit 9cfb841e38c7909da645b4c46b808fc94434951f. The source field is analysis_duration_sec: these are editable recording presets derived from analysis parameters, not a declaration of the full experimental protocol. Test mode uses seconds. Custom allows other durations. Each NPR recording is independent; an optional note can identify its day or phase. Presets are bundled for offline use and do not silently update from GitHub.

## Metadata and records

Mouse ID, cage, sex, group, weight, operator and notes are all optional. Record IDs use local date/time and a sequence, e.g. 20260922_143025_OFT_001_a7c9abcd. Filenames do not depend on free-text identifiers. There is no 12-mouse limit. Select a behavior or type an ID/date in the library filters.

Edit corrects a completed/interrupted/failed session and retains prior values and a reason. Delete hides its record and excludes it from the CSV log; it never deletes videos. Show deleted records and Restore undo it. JSON retains all records. Recording and finalizing sessions cannot be edited/deleted. Session IDs, acquisition times and QC are not changed by metadata edits.

Download log (CSV) opens in Excel; it includes recording and optional mouse information. Download all records (JSON) preserves events, revisions and quality metrics. User text beginning with spreadsheet formula markers is prefixed with an apostrophe in CSV exports for literal interpretation.

## Files and SongLabDLC

Output: project folder / data / TEST or EXPERIMENT / behavior / recording_ID /

- recording_ID.mkv: original H.264 capture with acquisition timestamps.
- recording_ID.mp4: stream-copy MP4, published only after decoded frame hashes/counts and relative timestamps match the original (1.1 ms tolerance).
- recording_ID_frame_timing.csv: zero-based frame index, relative timestamp, interval and decoded-frame MD5.
- session.json: full session metadata and quality report.
- recording_capture.log: FFmpeg/camera diagnostics.

Both video copies remain on disk (approximately twice the video storage). Processing can take several minutes for a long recording. Conversion failure leaves MKV; an unclean system shutdown may leave an unverified .pending.mp4; never use a pending file as a verified output. Early stop remains interrupted even when its retained video can be converted.

Run DLC on the intended MP4, then save the DLC CSV beside it. The MATLAB 0.2.0 companion matches basenames and can choose .mp4 when MKV is also present:

    SongLabDLC_behavior_analysis("OFT", "D:\Data\EXPERIMENT\OFT\recording_ID", ".mp4")

For batch analysis, copy chosen videos and matching DLC CSVs into one analysis folder per behavior. The analysis tool does not recursively gather nested recording folders. Keep each original session folder for archival metadata; do not mix behaviors or duplicate tracking outputs.

## Manual alignment

Start calibration, freeze a frame, select four rectangle corners clockwise from top-left (or drag ellipse bounds), then Track selected arena. Adjust the camera by hand while KLT features and a robust affine transform update the outline. Rotation, opposite-edge asymmetry, center offset and circle axis ratio provide visual guidance. Confirm alignment only after Alignment OK. This adapts the v2 workflow using OpenCV; it is not a bit-for-bit MATLAB tracker port or a physical calibration certificate.

Calibration targets 25 fps and displays measured camera delivery, tracking and browser refresh rates separately. It may run slower on your camera/computer. Recording automatically stops tracking and uses the lower-load preview. No overlays are burned into videos. Each behavior can save an outline in the project. Loaded profiles always require verification (including same-day reloads); freeze the current image, load the outline, check geometry and track again.

## Projects and idle processing

New Project creates a uniquely named folder containing BehaviorHub.project.json and data/. Save Project also saves the current recording settings. Records, events, corrections, queue changes and calibration selections autosave. Open Project restores records/settings and pending work. One project can span behaviors and days. Only one service may open a project at a time. Close it before copying the **whole project folder**; video references are relative. A project file alone does not contain videos. Recent projects are local shortcuts; after moving computers use Open Project at the new location.

Legacy Library preserves original records. Import legacy records (copy) explicitly copies available non-deleted old recordings into the current project; it never moves originals. Imported copies consume additional disk space.

MP4 remux and full frame QC run one job at a time, at reduced process priority, after the camera is closed and there has been no user activity for 60 seconds (adjustable 5-3600 seconds). Process pending videos starts work immediately. Pause, Resume and Retry failed manage the queue. Opening preview/calibration or starting recording cancels the background subprocess and returns its job to the queue; work restarts later rather than resuming at a byte offset. MKV is retained. Failed conversion does not invalidate or delete the source recording. The service must remain open to process videos; Exit application preserves pending work for the next project load.


## Timing quality

Preview produces thumbnails at up to 8 fps and is separate from the requested recording FPS; both outputs still share camera and CPU resources. A smooth preview does not prove capture quality. The application never silently duplicates frames to manufacture 25 fps.

Quality reports include decoded frames versus expected, observed FPS from timestamp span, maximum gap, intervals over 1.5 times the nominal interval and estimated missing frames from gaps. These are not hardware loss counters. Count deviations above 0.5%, long gaps, non-increasing timestamps or processing failures trigger REVIEW. A timing pass does not certify focus, visibility, exposure, lighting or scientific analysis suitability.

Existing SongLabDLC time-based metrics assume uniform nominal/average FPS. Container compatibility does not correct variable frame-rate timing. Validate the actual recording computer before experiments.

## Source and rebuilding

Source/ contains the Python service, page, preset source and tests. Python 3.11; imageio-ffmpeg 0.6.0; PyInstaller 6.22.3; NumPy 1.26.4; OpenCV headless 4.11.0.86. Put ffmpeg.exe beside app.py (or install imageio-ffmpeg) for source use. Build with the supplied BehaviorHub.spec. For DOM tests, npm ci installs the pinned jsdom dependency. The test_package.py script expects dist/1.1.0/BehaviorHub/BehaviorHub.exe. UI tests emulate a browser and do not replace a real-browser acceptance check.

## Upgrade from SongScope

Behavior Hub is the new product name. Finish recording and finalization, then click Exit application in SongScope before starting BehaviorHub.exe. Extract the new package into a separate folder. Existing records are reused from %LOCALAPPDATA%\SongScope; that internal directory name is retained for compatibility. The service port remains 43831. Opening the new EXE while SongScope is still running reconnects to that old service, so exit it first to see the new name. Video files and recording IDs are unchanged.

## Run or build from this repository (Windows)

The repository contains source code, not compiled binaries or experimental recordings. From this directory:

```powershell
python -m pip install -r requirements.txt
python app.py
```

To create `dist/1.1.0/BehaviorHub/BehaviorHub.exe` and the shareable ZIP:

```powershell
python build_windows.py
```

The build script copies the FFmpeg executable supplied by the pinned imageio-ffmpeg package, runs PyInstaller with BehaviorHub.spec, and packages the executable, dependencies, documentation, licenses and source. Run commands from this directory. No live experiment or camera is needed for the automated tests:

```powershell
python -m unittest test_behaviorhub test_projects -v
npm ci --ignore-scripts
python test_package.py
```

The packaged test uses a generated video source and isolated temporary records. Testing does not establish hardware capture quality; see VALIDATION.md. Runtime needs no MATLAB or DLC installation. Those remain separate analysis steps.
