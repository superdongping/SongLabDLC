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
