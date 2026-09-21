# Multi-format validation (2026-09-21)

Environment: Windows x64, MATLAB 24.1.0.2603908 (R2024a) Update 3, FFmpeg 7.1.

Passed test_video_formats with one-second 160x120/25fps generated fixtures:
- H.264 MP4, MKV and MOV; MJPEG AVI; exactly 25 frames decoded for each.
- Forced native-read failure path: stream-copy MP4, full decoded hashes/count and relative timestamp verification, readable by VideoReader.
- Same-basename containers raise AmbiguousVideo unless an extension is selected.
- Duplicate DLC tracking files raise DuplicateCSV; per-frame timing CSVs are excluded.
- MATLAB parsing checks for helper and assay files, plus main entry point.

Not tested: R2024b or R2025 releases, every codec, full real-data assay output equivalence, interactive ROI selection and the format-choice dialog in a real user session. Existing scientific algorithms and uniform-FPS assumptions remain unchanged.

The change is limited to the unified SongLabDLC_behavior_analysis workflow and its helpers/assays. Archived standalone scripts and fiber-photometry workflows are not migrated.

Run tests with addpath('tests'); test_video_formats(fixtureDir). The directory must initially contain sample.mp4, sample.mkv, sample.mov and sample.avi with 25 frames each, and no DLC test CSVs. The test writes fixture CSV files into that directory. Use a fresh disposable directory.
