# SongLabDLC

SongLabDLC contains MATLAB tools for behavior analysis after DeepLabCut (DLC), with a standardized pipeline for common mouse behavioral assays.

The main behavior-analysis workflow is now a single MATLAB entry point:

```matlab
SongLabDLC_behavior_analysis
```

This is the recommended way to run behavior analysis.

## Supported Behavior Tests

The standardized pipeline supports:

- Open Field Test (OFT)
- Novel Object Recognition / Novelty Preference Recognition (NPR)
- Zero Maze
- Y-maze
- Forced Swim Test (FST)
- Tail Suspension Test (TST)

## What You Need Before Running

You need:

- MATLAB
- DLC output CSV files
- Matching MP4 behavior videos
- The CSV and MP4 files for one assay placed in the same folder

For example:

```text
My_OFT_Data/
  Mouse1.mp4
  Mouse1DLC_..._filtered.csv
  Mouse2.mp4
  Mouse2DLC_..._filtered.csv
```

The pipeline tries to match CSV and MP4 files by file name. It is best if the video name appears at the beginning of the DLC CSV name.

Good example:

```text
2026-04-29_13-32-06.mp4
2026-04-29_13-32-06DLC_Resnet50_..._filtered.csv
```

Avoid mixing different behavioral assays in the same input folder. Run OFT, NPR, Zero Maze, Y-maze, FST, and TST in separate folders.

## How To Run The Pipeline

### Step 1: Open MATLAB

Open MATLAB.

### Step 2: Set The MATLAB Current Folder

Set the MATLAB Current Folder to the SongLabDLC repository folder.

Example:

```text
C:\Users\pingdong\Desktop\Antigravity_test\SongLabDLC
```

You should see `SongLabDLC_behavior_analysis.m` in the file list.

### Step 3: Run The Main Script

In the MATLAB Command Window, type:

```matlab
SongLabDLC_behavior_analysis
```

Then press Enter.

### Step 4: Select The Assay

A window will ask which assay you want to analyze.

Choose one:

- `OFT`
- `NPR`
- `Zero Maze`
- `Y-maze`
- `FST`
- `TST`

### Step 5: Select The Data Folder

Choose the folder containing the matching DLC CSV files and MP4 videos.

The pipeline will show a table of matched files in MATLAB. Check that each CSV is paired with the correct video.

### Step 6: Follow The On-Screen ROI Or Scaling Instructions

The pipeline will show the first video frame and ask you to draw or click assay-specific regions.

## Assay-Specific Instructions

### Open Field Test (OFT)

The pipeline will ask you to:

1. Draw the outer region for analysis.
2. Draw the open field box.
3. Draw a scale line representing `0.45 m`.

Outputs include:

- total distance traveled
- time in center
- center time ratio
- center distance
- total speed
- center speed
- trajectory QC image
- occupancy heatmap in seconds

The OFT heatmap uses cm on the X/Y axes and seconds on the colorbar. Heatmaps from the same batch use a shared color scale for easier comparison across mice.

### Novel Object Recognition / Novelty Preference Recognition (NPR)

The pipeline will ask you to:

1. Draw the open field box.
2. Draw a scale line representing `0.45 m`.
3. Click the center of object ROI 1.
4. Click the center of object ROI 2.

The object ROI radius is set to `0.04 m`.

Outputs include:

- object 1 exploration time
- object 1 event count
- object 2 exploration time
- object 2 event count
- discrimination ratios
- warning notes for low or zero exploration
- event timestamps
- trajectory and ROI QC image

NPR exploration uses bout-based event detection. Brief gaps are merged, and very short bouts are removed.

### Zero Maze

The pipeline will ask you to:

1. Draw the first closed-arm ROI polygon.
2. Draw the second closed-arm ROI polygon.

Outputs include:

- total distance in pixels
- time in open arms
- ratio of time in open arms
- time in closed arms
- trajectory QC image
- heatmap

### Y-maze

The pipeline will ask you to:

1. Draw the A-arm polygon.
2. Draw the B-arm polygon.
3. Draw the C-arm polygon.

Outputs include:

- arm-entry sequence
- spontaneous alternation percentage
- trajectory QC image

### Forced Swim Test (FST)

The pipeline will ask you to:

1. Draw a scale line for the beaker height, `24.25 cm`.

Outputs include:

- total immobilized time from `0-300 s`
- latency to first immobility
- total speed plots
- speed heatmap
- weighted immobility score plot

FST uses the original FST weighting scheme:

- Head: 5%
- MainBody: 5%
- Butt: 5%
- TailTip: 5%
- LeftForehand: 10%
- RightForehand: 10%
- LeftHindpaw: 30%
- RightHindpaw: 30%

### Tail Suspension Test (TST)

The pipeline will ask you to:

1. Draw a scale line for the beam width, `52 cm`.

Outputs include:

- total immobilized time from `0-360 s`
- latency to first immobility
- total speed plots
- speed heatmap
- weighted immobility score plot

TST uses the original TST weighting scheme:

- Head: 30%
- MainBody: 20%
- Butt: 10%
- TailTip: 0%
- LeftForehand: 10%
- RightForehand: 10%
- LeftHindpaw: 10%
- RightHindpaw: 10%

## Output Folder

Each analysis creates a folder inside your selected data folder:

```text
SongLabDLC_behavior_results/
```

Inside that folder, each run gets a timestamped subfolder.

Example:

```text
SongLabDLC_behavior_results/
  OFT_20260602_150748/
    summary_Open_field.xlsx
    run_metadata.json
    Mouse1_FirstFrame_with_Trajectory.tif
    Mouse1_SmoothedOccupancyHeatmap_sec.tif
```

## What To Check After Each Run

Always check the QC images before trusting the summary table.

Check:

- The CSV and MP4 files were matched correctly.
- The trajectory is inside the correct arena.
- The ROI or arm polygons are drawn correctly.
- The scale line is correct.
- The heatmap orientation matches the video.
- The summary Excel file has reasonable values.
- `run_metadata.json` records the run settings.

## Common Problems

### The wrong CSV is matched with the wrong MP4

Rename the files so the video name and DLC CSV name start with the same text.

Example:

```text
Mouse1.mp4
Mouse1DLC_..._filtered.csv
```

### MATLAB says a body part was not found

The DLC CSV body-part name may be different from what the pipeline expects.

Expected common names:

For OFT, NPR, Zero Maze, and Y-maze:

- Nose
- Head
- MainBody
- Butt
- TailTip

For FST and TST:

- Head
- MainBody
- Butt
- TailTip
- LeftForehand
- RightForehand
- LeftHindpaw
- RightHindpaw

### The heatmap looks too dim or too saturated

The OFT heatmap uses a shared batch color scale. This makes mice comparable within a run. If a cohort has extreme corner-staying behavior, the heatmap scale may need adjustment in:

```matlab
helpers/get_default_behavior_options.m
```

Look for:

```matlab
options.heatmap_color_limit_percentile = 99.5;
```

Lower values make heatmaps brighter but more saturated. Higher values make heatmaps less saturated but dimmer.

## Repository Organization

Recommended behavior-analysis files:

```text
SongLabDLC_behavior_analysis.m
assays/
helpers/
```

Older standalone behavior scripts are archived in:

```text
archived_behavior_scripts/
```

These archived scripts are kept as historical references and backups.

Other tools remain in:

```text
Fiber photometry V4_3/
Multi-fiber photometry V6_9/
Fiber_photometry_video_sync_tag_app_V1/
arena-live-alignment/
```

## Major Update Log

### June 2, 2026

- Added the standardized `SongLabDLC_behavior_analysis.m` entry point.
- Combined OFT, NPR, Zero Maze, Y-maze, FST, and TST into one user-facing workflow.
- Added shared helpers for DLC CSV loading, CSV/MP4 matching, video metadata, summary output, and run metadata.
- Added `run_metadata.json` output for reproducibility.
- Added batch-aware OFT occupancy heatmaps with cm axes and seconds-based colorbar.
- Corrected OFT heatmap orientation to match the video trajectory image.
- Tuned OFT heatmap color scaling for cohort visualization.
- Preserved original FST and TST assay-specific scaling, time windows, speed thresholds, and immobility weights.
- Moved older standalone behavior scripts into `archived_behavior_scripts/`.
- Moved the standardized pipeline entry point, `assays/`, and `helpers/` to the repository root.

### Previous Script Updates

- `NPR_post_DLC_V5_1_batch.m`: Added bout-based NPR exploration detection, merged short gaps between exploration frames, removed very short bouts, added discrimination ratios, and added warning notes for low or zero exploration.
- `NPR_post_DLC_V5_batch.m`: Updated NPR ROI definition by using a `0.45 m` scale line, user-clicked object centers, and a predefined object ROI radius of `0.04 m`.
- `Open_field_post_DLC_V3_batch.m`: Previous latest standalone OFT batch-processing script.
- `Zero_Maze_post_DLC_V2_batch.m`: Previous latest standalone Zero Maze script.
- `Y_Maze_post_DLC_V1_5.m`: Previous latest standalone Y-maze script for arm sequence and alternation percentage.
- `FST_DLC_V2.m`: Previous latest standalone FST script.
- `TST_DLC_V2.m`: Previous latest standalone TST script.

## DeepLabCut Resources

DeepLabCut documentation:

https://deeplabcut.github.io/DeepLabCut/docs/installation.html

DeepLabCut YouTube channel:

https://www.youtube.com/@deeplabcut7702

## Contact

For questions, contact superdongping@gmail.com or pingdong@unc.edu.
