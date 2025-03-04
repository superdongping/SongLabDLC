SongLabDLC
This repository contains MATLAB scripts designed for batch processing behavioral analysis data, specifically for use after extracting mouse body parts with DeepLabCut (DLC). 

If you have any questions about how to install and use Deeplabcut (DLC), find more information here:
https://deeplabcut.github.io/DeepLabCut/docs/installation.html

https://www.youtube.com/@deeplabcut7702

The repository focuses on processing data from various behavioral assays, including:

Open Field Test (OFT)
Novelty Preference Recognition (NPR)
Zero Maze
Forced Swim Test (FST)
Tail Suspension Test (TST)
DLC Extracted Body Parts
For different behavioral paradigms, the following mouse body parts are extracted using DLC:

OFT, NPR, Zero Maze:

Nose
Head
Main body
Butt
Tail tip

FST and TST:

Head
Main body
Butt
Tail tip
Left forehand
Right forehand
Left hindpaw
Right hindpaw

Instructions for Use

Use the Latest Script Versions:
For each behavioral assay, always use the most recent version of the script to ensure accuracy and updated features. For example:

For NPR, use NPR_post_DLC_V4_batch.m instead of NPR_post_DLC_V2_batch.m.

Similarly, for other behavioral assays, opt for the script with the highest version number.

Keep Files Organized:
Always keep the video files and their corresponding DLC-analyzed CSV files in the same folder as the MATLAB script. This ensures smooth processing without file path errors.

File Descriptions

DLC_Manual_Compare_NPR.m: Compares manual annotations with DLC outputs for NPR data.

FST_DLC_V2.m: Processes data for Forced Swim Test (FST).

NPR_post_DLC_V4_batch.m: Latest version for batch processing NPR data.

NPR_post_DLC_V5_batch.m changes the way the user to define the two ROIs for the objects. It asks the user to draw a line of one of the open field boarder as 0.45m. It changed to let the user click to specify the centers of the two ROIs. For the ROI radius, it is pre-define as 0.04 m.

Open_field_post_DLC_V3_batch.m: Latest version for batch processing Open Field Test data.

Zero_Maze_post_DLC_V2_batch.m: Processes Zero Maze data.

fiber_photometry_behavior_analysis.m: Synchronizes fiber photometry data with videos from behavioral events.

summarize_fear_data.m: Summarizes fear-related behavioral data.

Y_Maze_post_DLC_V1_5.m to process the Y maze data after DLC. The user need to draw 3 arms, arm A, arm B and arm C. Then the MATLAB code will quantify the Sequence and Alternation_Percentage from the Y-maze test.

If you have any questions, please get in touch with me at superdongping@gmail.com or pingdong@unc.edu
