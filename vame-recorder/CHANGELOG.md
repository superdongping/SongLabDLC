# 1.3.0 - 2026-09-25

- Persistent, cancellable idle queue for each baseline/post MP4 and full timing QC.
- Recording and injection workflow no longer wait for conversion; new capture/preview preempts processing.
- Automatic idle waits for completion/end of the current two-phase session, closed camera and 60 seconds of no app activity (configurable).
- Setup preview up to 25 fps; recording preview up to 8 fps, both width 640. Original recording FPS and timestamp policy unchanged.
- Queue pause/resume/retry/manual controls; MKV retained, restart recovery and existing-record compatibility.

# Version 1.2.0 (2026-09-21)

- Default resolution is now 1280 x 720 in both UI and service fallback.
- Edit completed/interrupted/failed session metadata with a correction reason and retained prior values.
- Delete records from the normal list and Excel without deleting video files; show deleted records and restore them.
- Active and waiting sessions are protected. Edits do not alter other sessions, roster data, capture times or QC.

# Version 1.1.0 (2026-09-21)

- Retain original MKV and automatically prepare MP4 after each phase without re-encoding.
- Verify decoded frame hashes/count and relative presentation timestamps before publishing the MP4 filename.
- Export per-frame timing CSV; report measured FPS, maximum gap, long intervals, frame deficit and estimated missing frames.
- Keep the UI responsive during processing; do not permit a new capture or shutdown during finalization.
- Account for both retained video copies in the initial storage budget.
- Keep acquisition status separate from timing quality and MP4 verification. Existing recordings remain untouched.

# Changes

## 1.0.1

- All application labels, notifications, errors and usage instructions are in English. Existing user-entered data is preserved verbatim.
- Reopening the executable connects to the existing VAME service instead of attempting to start a second server and failing on its data-directory lock.
- Windows sockets now use exclusive address binding. The data-directory lock is acquired before the service binds its port or reads the registry.
- Simultaneous launches wait for the first instance to become ready. A service using the same data directory on another port can also be reopened.
- Port conflicts with unrelated applications and other startup problems produce a concise English message instead of an unhandled-exception traceback.

## 1.0.0

- Initial local recording service, two-phase acquisition workflow, dose calculator and printable Excel logs.
