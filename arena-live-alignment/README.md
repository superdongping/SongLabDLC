# Arena Live Alignment Unified

A MATLAB tool for **live camera alignment** for behavioral arenas using a USB webcam.  
It supports both:

- **Rectangular arenas** (for example, open field box)
- **Circular arenas** (for example, zero maze)

The program helps users adjust the overhead camera so that the arena is as parallel and centered as possible before recording.

## Features

- Automatically detects connected webcams
- Lets the user select which webcam to use
- Supports two arena types:
  - **Rectangular arena**
  - **Circular arena**
- Starts with **live preview**
- Freeze frame for calibration
- Dynamic live feedback during camera adjustment
- Recalibration without restarting the program

## Requirements

- MATLAB
- MATLAB Support Package for USB Webcams
- Computer Vision Toolbox

## File

- `arena_live_alignment_unified.m` — main program

## How to Run

Save the script in your MATLAB working directory and run:

```matlab
arena_live_alignment_unified
