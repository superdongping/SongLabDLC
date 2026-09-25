# VAME Recorder 1.3.0 - 2026-09-25

Synthetic regression coverage: 16 distinct service/media tests covering original two-phase gates, dose arithmetic, metadata/export protections, interruption/recovery, deferred phase processing, no automatic work during injection wait, next-animal preemption, persistent paused queue, automatic idle processing, failed conversion retry, source preservation, true FFmpeg subprocess cancellation, preview delivery above the old rate and recording preview target 8 fps. Packaged HTTP/DOM/relaunch tests verify the final executable separately.

No physical camera was opened for these tests. No new two-hour stability, real-browser visual or native picker validation is claimed. Higher preview targets do not prove camera acquisition FPS; earlier physical-camera QC limitations remain. Saved video timing must be checked after idle processing completes. Validation files use isolated state, not the active experiment registry.

## Historical validation

# Version 1.2.0 validation (2026-09-21)

- 12 automated tests passed, including 720p default, invalid-edit rejection, protection of active/waiting sessions, independent session corrections, recalculated volume, retained audit history, delete/restore persistence and unchanged original video bytes.
- Additional export assertions passed: deleted experiment omitted from printable Excel; restored experiment included again.
- Packaged EXE DOM-emulation flow passed: default 1280x720, two phases, MP4/QC display, edit dialog and saving, delete from list, show deleted and restore. Dialog mechanics were emulated; real-browser visual acceptance is still outstanding.
- Packaged HTTP flow and full relaunch regression passed. No real camera was opened for these synthetic tests, and no existing experimental records were edited.
- A 720p default does not establish 25fps acquisition quality. Actual target-computer camera timing and clarity qualification remain necessary.

# Version 1.1.0 validation ? 2026-09-21

- 11 automated tests passed: existing two-phase workflow, recovery/interruption and new conversion/QC tests. Normal phases produce both retained MKV and verified MP4. A deliberately omitted frame preserves its timestamp gap and triggers REVIEW. Corrupt input, conversion failure and existing destination do not falsely report verified output or overwrite originals.
- Packaged Windows EXE HTTP tests passed: authentication, recording, automatic stop, both phases, refresh/state, saline, Excel export and shutdown.
- Packaged EXE relaunch tests passed: repeated opening, alternate port, reopening during recording, clean restart, concurrent cold launch, unrelated occupied port and preserved user data.
- DOM emulation passed: English forms, dose, folder entry, test durations, both phases, reload, controls, measured 25 FPS and verified MP4 display. This is NOT a real-browser visual or native folder-picker acceptance test.
- Existing real C920 10-second sample was COPIED to validation/mp4_1.1.0 for processing. Results: 227/250 frames; measured 22.690763 FPS; maximum gap 120 ms; 20 intervals above 60 ms; estimated 23 missing frames from gaps. Correctly marked REVIEW. MP4 has identical decoded frame content/count and zero measured relative timestamp difference; original recording unchanged.
- The earlier two-hour camera run was still active when checked (~54 min, approximately 22.4 fps). This release does not claim completed full-duration validation or resolve camera acquisition quality.
- Full-duration processing performance, actual browser/folder-picker interaction and recording-computer exposure/lighting/focus/USB/frame-rate qualification remain required. TIMING_CHECKS_PASSED is a software timing screen, not a scientific acceptance certificate.

## Earlier validation

# Validation — VAME Recorder 1.0.1

Date: 2026-09-21. Host: Windows x64; Python 3.11.4; FFmpeg 7.1; Logitech HD Pro Webcam C920.

## Version 1.0.1 regression checks

- Repeated launches of the packaged executable exit cleanly and leave the original service alive, with the registry hash unchanged.
- Relaunch during a simulated recording preserves the session and all 150 expected frames in the six-second test.
- Exit and restart preserves saved mouse information, including Unicode user notes.
- Two simultaneous cold launches produce one service and one clean handoff; a shared data directory using another port also hands off correctly.
- An unrelated application occupying the port produces a controlled failure and does not create or overwrite a mouse registry.
- Six recorder unit tests and the packaged HTTP workflow pass. English DOM workflow and application-text checks cover forms, prompts, status messages and the printable template. DOM emulation is not a browser layout inspection.
- New application-owned strings and documentation are English. Historical user records are preserved verbatim, including their original language.

## Completed

- Unit tests: dose conversion at 20/25/30 g; invalid/nonfinite weights; stage ordering; duplicate-start/injection rejection; timer-driven two-stage recording; explicit interruption; restart recovery; literal text IDs in Excel.
- Packaged EXE HTTP test: starts without invoking system Python; loopback token checks; mouse and folder saving; saline calculation; two timed stages; page reload does not stop capture; Excel contains 13 sheets with cached dose values and print layout; orderly service exit. Video for this test is generated, not a camera.
- DOM emulation: mouse form, 0.3125 mL at 25 g, folder entry, seconds/minutes test toggle, saline injection prompt, two-stage controls and refresh recovery. No JavaScript DOM errors detected. This is not a real-browser rendering test.
- Spreadsheet: artifact-tool formula/link recalculation at 20/25/30 g; no formula errors; all 13 sheets visually inspected; Letter print areas and one-page-per-mouse page settings added and checked structurally. No physical print or native Excel UI test performed.
- Real C920, MJPEG, 1920x1080, requested 25 fps: baseline and post test files auto-saved and decoded without decoder failure. 10-second baseline: 227 frames, PTS 0–9.96 s; 15-second post: 344 frames, PTS 0–14.96 s. Maximum observed timestamp gap 120 ms. Files preserve these gaps rather than inserting duplicate frames.

## Not yet passed

- **25 fps frame completeness is not passed on the present physical setup.** Counts below 250/375 expected frames are flagged for review. Camera/USB/exposure/lighting and target-computer acquisition must be checked before claiming paper-matched constant-rate capture.
- Real-browser visual/button testing and the native directory chooser: browser-control tools returned no available browsers/apps. DOM and HTTP checks do not replace this final interaction check.
- Actual mouse roster is not supplied. The printable workbook deliberately contains 12 blank entries; CAMERA_TEST/DOM_TEST are synthetic validation identities only.

## Two-hour real-camera run

Started approximately 09:38 America/New_York on 2026-09-21, using a 10-second test baseline followed by a 7200-second post-return phase. It is isolated under `validation/TEST/20260921_093806_bf7fc8ac`. The injection event is simulated; no animal or actual injection is involved.

**Pending at initial packaging.** The recorder will stop automatically, then the validation script decodes both videos and writes `validation/camera_20260921_093806.json`. A separate finite watcher writes `LONG_RUN_RESULT.md` and updates the same local task-board task when the result is available. A successful long run still does not clear the frame-completeness warning above.

The target Windows PC, its storage and its actual camera setup still require their own short workflow and full-duration validation.
