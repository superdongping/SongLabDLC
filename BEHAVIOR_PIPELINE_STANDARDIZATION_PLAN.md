# Behavior Pipeline Standardization Plan

## Goal

Standardize the MATLAB post-DeepLabCut behavior analysis workflow so that the analysis is more reliable, reproducible, and easier to use.

The first priority is to consolidate the behavior-analysis scripts for:

- Open Field Test (OFT)
- Novel Object Recognition / Novelty Preference Recognition (NPR)
- Zero Maze
- Y-maze
- Forced Swim Test (FST)
- Tail Suspension Test (TST)

Fiber photometry scripts should remain unchanged during this phase.

## Current Problem

The repository currently contains multiple script versions for the same behavioral assay. For example, NPR has versions from `NPR_post_DLC_V2_batch.m` to `NPR_post_DLC_V5_1_batch.m`. This makes it hard for new users to know which script to run and increases the risk of inconsistent analysis.

Several scripts also repeat the same workflow steps:

- Find DLC CSV files
- Find corresponding MP4 videos
- Load DLC coordinates
- Read video frame rate
- Ask the user to define ROIs or scaling
- Apply likelihood filtering
- Run assay-specific analysis
- Save summary tables and QC figures

These common steps should be shared across assays.

## Proposed Direction

Create one main behavior-analysis entry point, for example:

```matlab
SongLabDLC_behavior_analysis.m
```

The user would run this one file, then select:

1. Behavioral assay type
2. Folder containing DLC CSV and MP4 files
3. Analysis options, if needed

The main script would then call assay-specific analysis modules.

## Recommended Repository Structure

```text
SongLabDLC/
  behavior_pipeline/
    SongLabDLC_behavior_analysis.m
    assays/
      analyze_oft.m
      analyze_npr.m
      analyze_zero_maze.m
      analyze_y_maze.m
      analyze_fst.m
      analyze_tst.m
    helpers/
      match_csv_video_files.m
      load_dlc_csv.m
      get_video_info.m
      save_run_metadata.m
      save_qc_figure.m
      validate_dlc_bodyparts.m
  archive_old_behavior_scripts/
  Fiber photometry V4_3/
  Multi-fiber photometry V6_9/
  README.md
```

This structure separates the new standardized behavior pipeline from older scripts and from fiber photometry code.

## Latest Behavior Scripts To Use As Starting Points

- OFT: `Open_field_post_DLC_V3_batch.m`
- NPR: `NPR_post_DLC_V5_1_batch.m`
- Zero Maze: `Zero_Maze_post_DLC_V2_batch.m`
- Y-maze: `Y_Maze_post_DLC_V1_5.m`
- FST: `FST_DLC_V2.m`
- TST: `TST_DLC_V2.m`

Older versions should not be deleted at first. They can be moved to an archive folder after the new pipeline is validated.

## Shared Components To Build First

### 1. File Matching

Create a shared helper to match DLC CSV files with MP4 videos.

This should be more reliable than assuming:

```matlab
csv_files(i) matches video_files(i)
```

The helper should:

- List all CSV files
- List all MP4 files
- Match by shared base name when possible
- Show a confirmation table
- Warn if files cannot be matched

### 2. DLC CSV Loader

Create a shared DLC CSV loader that can read standard DeepLabCut CSV files with multi-row headers:

- scorer
- bodyparts
- coords

The loader should return coordinates by body part name rather than fixed column number. This reduces errors when body-part order changes.

Example desired access pattern:

```matlab
dlc.Nose.x
dlc.Nose.y
dlc.Nose.likelihood
dlc.MainBody.x
dlc.MainBody.y
dlc.MainBody.likelihood
```

### 3. Video Metadata

Create a shared helper to read video information:

- frame rate
- duration
- width
- height
- number of frames, if available

### 4. Run Metadata Output

Every analysis run should save metadata, such as:

- assay type
- date and time
- script or pipeline version
- CSV file names
- MP4 file names
- CSV/video matched pairs
- frame rate
- likelihood threshold
- analysis duration
- ROI coordinates
- scale factor
- body parts used
- output files generated

Possible output formats:

- `run_metadata.json`
- `run_metadata.xlsx`
- or both

### 5. QC Outputs

Keep the current strength of the repo: visual QC.

Each assay should continue saving:

- trajectory overlay
- ROI overlay
- heatmap, when relevant
- summary Excel file

## GitHub Workflow

Use a branch for this work.

Recommended branch name:

```text
standardize-behavior-analysis
```

Workflow:

1. Create the branch from `main`.
2. Build and test the standardized behavior pipeline on the branch.
3. Validate with demo data and at least one real dataset per assay.
4. Open a pull request into `main`.
5. Merge only after validation.

This keeps `main` stable while the new workflow is being developed.

## Suggested Implementation Stages

### Stage 1: Planning and Branch Setup

- Create this planning document.
- Create a new Git branch.
- Decide the exact folder structure.

### Stage 2: Shared Helpers

- Build file matching helper.
- Build DLC CSV loader.
- Build video metadata helper.
- Build metadata-saving helper.

### Stage 3: Unified Entry Script

- Create the main behavior pipeline entry script.
- Add assay selection.
- Add folder selection.
- Add file matching preview.

### Stage 4: Assay Migration

Migrate assays one at a time:

1. OFT
2. NPR
3. Zero Maze
4. Y-maze
5. FST
6. TST

For each assay:

- Preserve existing output metrics.
- Preserve QC figures.
- Add metadata output.
- Compare new output against old script output.

### Stage 5: Cleanup

- Update `README.md`.
- Move old behavior scripts to an archive folder.
- Keep fiber photometry scripts unchanged.
- Add simple demo instructions.

## Validation Checklist

For each assay, confirm:

- CSV and MP4 files are correctly matched.
- DLC body parts are detected correctly.
- Frame rate is read correctly.
- ROI/scaling selections are saved.
- Summary Excel output matches expected values.
- QC images are generated.
- Metadata file is generated.
- Results are comparable to the existing latest script.

## Notes

This plan focuses only on behavior analysis after DeepLabCut. Fiber photometry code is intentionally excluded from this first standardization phase.
