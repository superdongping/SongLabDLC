# SongLabDLC

MATLAB tools for post-DeepLabCut behavioral analysis and related Song Lab workflows.

## Standard Behavior Analysis

For behavior analysis after DLC, run the unified entry point:

```matlab
SongLabDLC_behavior_analysis
```

This pipeline supports:

- Open Field Test (OFT)
- Novel Object Recognition / Novelty Preference Recognition (NPR)
- Zero Maze
- Y-maze
- Forced Swim Test (FST)
- Tail Suspension Test (TST)

The pipeline asks the user to:

1. Select the assay type
2. Select the folder containing DLC CSV files and MP4 videos
3. Define assay-specific ROIs or scaling lines

Each run saves a timestamped output folder containing:

- summary Excel file
- QC figures
- `run_metadata.json`

The shared analysis code is organized in:

```text
assays/
helpers/
```

## Archived Behavior Scripts

Older standalone behavior-analysis scripts have been moved to:

```text
archived_behavior_scripts/
```

These files are kept as historical references and backups. For routine analysis, use `SongLabDLC_behavior_analysis.m`.

## Other Tools

Fiber photometry and arena-alignment tools remain in their existing folders:

```text
Fiber photometry V4_3/
Multi-fiber photometry V6_9/
Fiber_photometry_video_sync_tag_app_V1/
arena-live-alignment/
```

## DeepLabCut Resources

DeepLabCut documentation:

https://deeplabcut.github.io/DeepLabCut/docs/installation.html

DeepLabCut YouTube channel:

https://www.youtube.com/@deeplabcut7702

For questions, contact superdongping@gmail.com or pingdong@unc.edu.
