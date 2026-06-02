# SongLabDLC Standardized Behavior Pipeline

This repository contains a standardized post-DeepLabCut behavior analysis workflow.

Run the unified entry point from MATLAB:

```matlab
SongLabDLC_behavior_analysis
```

The program will ask for:

1. Assay type
2. Folder containing DLC CSV and MP4 files
3. Assay-specific ROI or scaling information

Supported assays:

- OFT
- NPR
- Zero Maze
- Y-maze
- FST
- TST

The pipeline writes results to:

```text
SongLabDLC_behavior_results/
```

Each run creates a timestamped output folder containing:

- summary Excel file
- QC figures
- `run_metadata.json`

Older standalone behavior scripts are archived in `archived_behavior_scripts/` and kept as historical references.
