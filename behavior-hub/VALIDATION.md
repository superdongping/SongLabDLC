# Behavior Hub 1.1.1 alignment rollback - 2026-09-22

Only camera alignment was rolled back to the 1.0.1 static guides and temporary reference overlay. Portable projects, date-first unique filenames and recording-priority idle MP4/QC queue remain. Existing dynamic calibration profiles are preserved but inactive. No OpenCV/NumPy tracking runtime is bundled. The user reported successful testing and authorized GitHub publication on 2026-09-22.

Validation: 11 service tests cover recording/project/queue/management regressions. Packaged EXE HTTP + DOM test covers recording, background verification, reference capture/clear, static circle control, record edits/deletion/restoration, exports, duplicate launch and project reopen. Synthetic images and mocked native folder selection do not certify real camera performance or actual browser layout. The user accepted this build after testing. This feedback does not establish quantitative camera FPS, full-duration performance or MATLAB scientific regression coverage.

## Historical validation (previous versions only)

# Behavior Hub 1.1.0 validation - 2026-09-22

Passed on this Windows computer using synthetic video / isolated temporary projects:

- 12 Python tests: six behavior presets, optional metadata, timed capture, retained original MKV, verified remux and per-frame timing; project isolation and exclusive lock; relative paths after whole-folder copy; saved settings; consecutive animals without MP4 wait; durable queue restart; conversion failure and retry; automatic idle start and pause; background preemption before capture; actual FFmpeg subprocess cancellation; KLT translation tracking and loss invalidation; record corrections/deletion/restoration.
- Packaged EXE starts with bundled NumPy/OpenCV. Repeat launch reuses service. HTTP + DOM test covers project creation, six presets, default 720p, recording, queued then verified MP4, ellipse controls, record edit/delete/restore, English text, exports and project reopening after service restart.
- Additional packaged HTTP check starts calibration, freezes/selects an arena, then starts capture and verifies tracking is disabled.
- All capture tests used generated frames, not a physical camera. Native directory/file dialogs are mocked in the DOM test. No browser-control surface was available for visual testing.

Remaining acceptance: real browser layout and native pickers; target webcam exposure/lighting/focus and actual delivered FPS; manual tracking usability, 25 fps calibration performance and recovery under camera movement; full experiment-duration capture and processing responsiveness on target hardware. Earlier C920 short-video timing did not establish 25 fps. No new physical-camera or two-hour result is claimed.

MP4 verification proves frame hash/count and relative-timestamp equivalence to MKV, not absence of losses during camera acquisition. Background cancellation restarts the job later; it is not byte-level resume. Software guidance thresholds are not a physical camera calibration certificate. MATLAB R2024b/R2025 and scientific six-assay regression remain unverified; MATLAB time metrics still use nominal/average FPS.

## Previous release validation (historical)

# Behavior Hub 1.0.1 rename validation

The product was renamed from SongScope to Behavior Hub. The packaged EXE passed the existing HTTP/DOM recording, management, export and restart flow. An isolated upgrade test launched SongScope, recorded a synthetic video, reopened BehaviorHub.exe while the old service was running, then shut down and opened the new executable. The existing session object and original video SHA-256 remained unchanged. The new page title was Behavior Hub. The legacy internal application identifier, request header, port and %LOCALAPPDATA%/SongScope directory are intentionally retained for compatibility.

# Validation — 2026-09-21

Passed:
- Four Python test groups covering six assays with generated frames, defaults imported from source, blank optional metadata, 720p default, automatic stop, retained MKV and verified MP4, unique IDs, invalid inputs, interrupted capture, concurrent-start protection, record edits, deletion/restoration, restart persistence, CSV literal-text protection and unchanged video bytes.
- Packaged EXE HTTP + DOM flow: all six duration presets, default 720p, empty metadata recording, measured 25fps synthetic stream, MP4 verification, manual guide control, record edit/delete/restore and English UI text.
- Packaged repeated launch, clean shutdown/restart and CSV/JSON exports.
- Separate MATLAB companion: Windows R2024a Update 3 read H.264 MP4/MKV/MOV and MJPEG AVI (25 frames each). Forced MP4 stream-copy fallback passed hash/count/timestamp comparison. Ambiguous same-name containers and duplicate tracking CSVs rejected; timing CSV excluded. MATLAB syntax checks passed.

Not established:
- Real-browser appearance, native folder picker and alignment usability. The available browser-control inventory was empty; DOM emulation uses a stub canvas and dialogs.
- A real-camera test of BehaviorHub 720p/25fps, long-duration finalization performance or freedom from dropped capture frames. No experimental camera session was started or interrupted to run these software tests.
- MATLAB R2024b/R2025a/R2025b tests or every codec within each recognized extension.
- Complete six-assay scientific-result regression with real DLC datasets and manually defined ROIs. Scientific algorithms were not changed; video discovery/reader entry points were updated.

Known limitation: existing MATLAB metrics use nominal/average FPS; recordings with irregular timing require analysis review. Manual guides are not the optical-flow tracking implementation in arena-live-alignment.
