# Behavior Hub 1.2.3 user guide and publication - 2026-09-30

Added a yellow button with bold red User guide text in the header, immediately before the service status. It opens a compact, seven-step English guide with a Close button and native Escape support. Camera/recording behavior is unchanged from 1.2.2. The packaged EXE HTTP/DOM suite passed, including opening/closing the guide, recording, MP4/QC, exports and restart. Focus interaction tests also passed. The user accepted the update and authorized GitHub publication after adding this guide. That acceptance is not a quantitative FPS, optical-sharpness, or long-duration certification.

# Behavior Hub 1.2.2 simplified focus UI - 2026-09-30

Removed manual capability-read/apply/save buttons, numeric entry, technical focus status text, and the 1:1 detail view. Selecting Manual obtains capabilities automatically; slider input is debounced for 300 ms and applied to preview automatically. Successful latest settings are saved to the current project's device profile. Pending/in-flight changes temporarily disable recording startup; the backend still applies/verifies the recording connection and blocks failed verification. Camera capture transport is unchanged; applying a focus change still briefly restarts preview.

Focus DOM tests passed: controls removed, automatic capability discovery/apply/save, burst coalescing, newest slider value retained during a slow request, rejected changes not saved, recording lock and device-profile restoration. Packaged 1.2.2 EXE HTTP/DOM recording, idle MP4/QC, guide/reference, edit/delete/restore, export, duplicate-launch and restart tests also passed. No physical-camera test was run for this UI-only update, so the user's active camera session is not disturbed. Prior hardware results below apply to the unchanged capture implementation.

# Behavior Hub 1.2.1 icon update - 2026-09-30

User-supplied transparent PNG preserved in assets/behavior-hub.png and converted to a multi-resolution Windows ICO (16/24/32/48/64/128/256). The EXE and browser favicon use the same artwork. Focus and recording behavior are unchanged from 1.2.0. Previous package is retained. Verified all seven embedded EXE icon resources against the ICO images. Packaged EXE startup, version 1.2.1 and browser favicon HTTP response passed in isolated synthetic state; no physical-camera settings were changed.

# Behavior Hub 1.2.0 focus control - 2026-09-30 (local test release)

- 11 existing service/project/queue regression tests passed.
- 6 focus tests passed: malformed/out-of-range/mismatched readback rejection; no encoder/video on focus failure; device-scoped saved profiles and project reopen; recording control lock; duplicate camera names; lossless NUT relay with variable-frame timestamp/frame-hash equivalence.
- Physical HD Pro Webcam C920 smoke test: focus range 0-250, step 5, both modes supported. Three full-resolution preview -> manual-focus recording transitions passed; each short recording returned manual flag/value from the active capture filter. An out-of-range native request was rejected before creating video. Original autofocus mode restored after testing. Disposable test videos were removed.
- Each 2-second C920 smoke recording contained 48 frames, approximately 23.98 measured fps despite a requested 25 fps. Existing timing QC correctly returned REVIEW. All three MKV -> MP4 conversions verified frame/timestamp equivalence. This is NOT a 25-fps or loss-free acquisition claim.
- Packaged EXE HTTP/DOM test passed: normal recording, idle MP4 verification, guides/references, record edits/delete/restore, exports, duplicate launch, service restart and project reopen.
- Focus DOM test passed: device range/step, manual controls and outgoing request, saved profiles, invalidation after editing, controls locked during recording, and separate identity for same-name cameras.
- Packaged EXE + physical C920 API test passed: manual-focus recording, save, complete service restart, project reopen, saved manual profile applied and freshly verified, and invalid-focus rejection with no video file. Original autofocus mode restored.
- Live browser visual QA was unavailable: the browser-control inventory was empty. DOM tests cover interaction logic, not actual layout. No test processes were left running.

Pending user acceptance: optical sharpness at the experimental plane, longest normal experiment duration, three real-use transitions, USB unplug/replug recovery, and restart/settings retention in normal use. Device readback verifies settings only; it is not continuous monitoring or proof of sharpness. GitHub publication has not been authorized for this release.

## Historical validation (earlier versions only)

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
