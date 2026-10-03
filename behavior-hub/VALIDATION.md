# Behavior Hub 1.2.6-r2 publication - 2026-10-03

The user tested the local build, reported it working well and explicitly authorized GitHub publication. The tested EXE and runtime are unchanged; release preparation updates documentation and packaging only. Prior test results and their limitations remain below. User acceptance is not a quantitative latency, frame-rate or optical-sharpness certification.

# Local 1.2.6-r2 smooth framing update (2026-10-03)

- Zoom/pan/flip-only updates retain the same live camera/FFmpeg process and focus connection. Browser CSS applies the even-pixel crop and flips immediately, before the debounced persistence request. Full-field JPEGs continue during both setup and recording; the recorded MKV still uses the FFmpeg crop/scale/flip filter and MP4 is verified against it. Setup JPEG rate increased from 4 to 8 fps; recording preview remains 8 fps at reduced resolution.
- `python -m unittest test_video_transform test_behaviorhub test_projects test_focus test_video_output -v`: all 26 cases passed. New coverage keeps the same process/thread and nonempty advancing preview through 30 updates with no project; invalid transforms leave the process intact and a resolution change still restarts it. Browser crop and flip coordinates match encoder geometry for 80 combinations including fractional zoom and edge crops. Real FFmpeg quadrant recordings validate the saved transform, output size, MP4 verification and restored project settings.
- `node test_focus_ui.mjs`: passed, including focus failure handling and recording locks.
- Isolated source service tested in the in-app browser with a simulated camera: zoom 4x, pan right/down, horizontal flip and Reset; decoded 1280-pixel preview stayed visible with the expected CSS transform. No user's project or real camera was opened. This is not a quantitative real-webcam latency or frame-rate certification.
- Updated one-page PDF rendered and visually checked. Packaged EXE HTTP/DOM regression passed: immediate CSS update before debounce, reference transform parity, Reset/latest slider value/pan direction, recording locks, timed recording, MP4/QC, Auto_ID guide, export, project settings and restart. The first new DOM assertion rounded 1.300813 to 1.301 incorrectly; corrected its expected value to the exact source/crop ratio and reran successfully. Final ZIP inventory/hashes are verified separately at packaging.

# Behavior Hub 1.2.6 zoom and flip - 2026-10-02 (local test build)

Added digital zoom (1-4x), pan, horizontal/vertical flips and Reset. The same validated FFmpeg framing filter applies to both preview and the recorded MKV, preserving selected output dimensions; verified MP4 remains a stream copy of that transformed MKV. No uncropped copy is retained. Settings save with the project and session. UI changes debounce, clear references and restart an open preview; recording locks the controls. Pan directions account for flips. Changing camera/resolution resets UI framing.

The 24-test Python suite passed after updating a focus-test fixture to include the new default transform fields (the initial old fixture failed its strict live-settings comparison). Real FFmpeg tests use a static four-color source to verify preview/recording orientation for each flip combination and zoom/pan, output dimensions, MP4 frame/timestamp equivalence, invalid settings and restart persistence. Focus DOM and packaged EXE HTTP/DOM regression passed, including Reset, coalesced zoom updates, flipped pan direction, recording locks, transformed capture, MP4/QC, saved settings, export and reopen. PDF guide rendered and visually inspected. No physical camera was opened. Real-camera throughput, experimental framing and browser appearance remain for user acceptance. No GitHub upload is authorized for this build.

# Behavior Hub 1.2.5-r4 publication - 2026-10-02

The user tested the local 1.2.5-r4 build, reported it looked good, and explicitly authorized GitHub publication. The executable and application assets are unchanged from that tested build. Release preparation updates documentation and the PDF publication label only. Prior test results and quantitative camera limitations below remain applicable.

# Behavior Hub 1.2.5-r4 Auto_ID labels - 2026-10-02 (local test build)

Blank Mouse ID now produces Auto_ID01, Auto_ID02, etc. Existing project counters continue without resetting; old video filenames remain unchanged. Updated the naming hint, top User guide and one-page PDF to explain that Auto_ID is an automatically assigned mouse ID when the field is blank, with numbering saved per project across restarts.

Five targeted output tests and packaged EXE HTTP/DOM regression passed, including automatic-ID filenames, numbering across restart, custom behavior, MP4/QC, guide wording, exports and saved-project reload. The one-page PDF was rendered and visually checked. No physical-camera validation or GitHub publication was performed.

# Behavior Hub 1.2.5-r3 custom behavior and filenames - 2026-10-02 (local test build)

Custom now requires a behavioral test name. New MP4 filenames always use date/time, behavior name and Mouse ID. Blank IDs use a project-persisted counter (01, 02, 03), advanced only for blank-ID recording attempts. TEST prefix and collision suffixes remain. Names freeze at recording start. Existing sessions and filenames are unchanged. Custom names appear in the library, details and CSV log, while the internal CUSTOM category stays compatible with existing project folders.

22 Python tests, focus DOM tests and packaged EXE HTTP/DOM regression passed. Tests include blank custom-name rejection, custom-name sanitization, new filenames, CSV labels, shared output, persistent numbering across restart and named/unnamed recordings, project reopening with restored custom name, MP4/QC and export. In-app and one-page PDF guides updated; PDF rendered and visually checked. No physical camera was opened and no GitHub upload was performed. User acceptance is still required for actual workflow and browser layout.

# Behavior Hub 1.2.5-r2 UI revision - 2026-10-02 (local test build)

Moved the countdown into the recording controls row beside Stop early and save, with yellow fill and bold red text. It uses captured media time, displays Saving during finalization and Ready when idle. Moved Idle processing queue below Record library and replaced the stray question mark in the heading. View now opens a details dialog with an explicit Open video action for verified MP4s, including older session-folder outputs. Video paths come from the selected session, never a client-supplied filename. No historical files are relocated.

Targeted output/playback tests passed for shared and legacy destinations, blocked unverified/missing files and invalid paths. OS player launch was mocked to avoid opening applications during tests. Focus DOM and packaged EXE HTTP/DOM regression passed, including countdown placement/colors/font weight, live remaining-time display and idle reset, section order, View dialog, playback request, MP4 verification, export and project restart. Tests use isolated synthetic video, not the physical camera. Real Windows player association and final browser layout remain for user acceptance. Read-only inspection of the user's screenshot project showed its six recordings were made with 1.2.4, explaining why those MP4s use the old per-session layout.

# Behavior Hub 1.2.5 shared MP4 output - 2026-10-01 (local test build)

New recordings publish verified MP4 files into the project's MP4 folder. Users choose recording ID or date/time plus Mouse ID before recording. Names freeze at recording start, TEST names carry a prefix, and collisions receive numeric suffixes. Original MKV and timing reports stay in session folders. Existing sessions keep their original paths. No migration or GitHub upload was performed.

21 Python tests passed, including shared-folder output, sanitization, case-insensitive collisions, path escape rejection, fixed filenames after metadata edits, project relocation, conflict preservation, retry/restart, and verification of an already-published file. Existing capture/project/focus tests remain passing. The focus DOM suite and packaged EXE HTTP/DOM capture, MP4/QC, new filename UI, export and restart tests passed. The updated one-page PDF was rendered and visually inspected.

No physical camera session was opened. This is a local test candidate pending user validation of actual workflow and output files. Earlier C920 frame-rate and optical-quality limitations remain applicable. Keep 1.2.4 separately and exit it before starting 1.2.5.

# Behavior Hub 1.2.4 header icon - 2026-09-30

Added the existing transparent application artwork to the top-left header beside the Behavior Hub title, including a smaller layout on narrow screens. Original 1.2.3 release and PDF guide are retained. Recording and focus behavior are unchanged. Packaged EXE startup/version checks, the PNG route (byte-identical to the supplied artwork), and header DOM structure checks passed. No physical camera session was opened for this cosmetic update.

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
