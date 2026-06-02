clear all;   
close all;
clc

% Set the likelihood threshold (adjust as needed)
likelihood_threshold = 0.03; % Lowered from 0.05 to detect more events

% Search for CSV and MP4 files in the current folder
csv_files = dir('*.csv');
video_files = dir('*.mp4');

% Ensure there is a one-to-one correspondence between CSV and MP4 files
if length(csv_files) ~= length(video_files)
    error('The number of CSV files does not match the number of MP4 files.');
end

if isempty(csv_files) || isempty(video_files)
    error('No CSV or MP4 files found in the current directory.');
end

% Load the first video to define the open field and ROIs
first_videoFileName = video_files(1).name; 
fprintf('Loading first video: %s\n', first_videoFileName);

% Load the first video and get the first frame
videoObj = VideoReader(first_videoFileName);
firstFrame = readFrame(videoObj);

% Display the first frame and let the user draw a rectangle for the open field area
figure;
imshow(firstFrame);
title('Draw a rectangle to define the open field area');
openFieldRect = drawrectangle('Color','g'); % User defines the open field area

% Get the position and dimensions of the rectangle
openFieldPosition = openFieldRect.Position; % [x, y, width, height]

% -------------------------------------------------------------------------
% NEW ROI DEFINITION BASED ON A RULER AND MOUSE-SELECTED CENTERS
% -------------------------------------------------------------------------

% 1. Ask the user to draw a ruler to define a known length (0.45 m)
title('Draw a line representing a 0.45 m ruler');
rulerLine = drawline('Color','y'); % User draws a line
rulerPos = rulerLine.Position;    % [x1, y1; x2, y2]
rulerPixelLength = sqrt( (rulerPos(2,1)-rulerPos(1,1))^2 + (rulerPos(2,2)-rulerPos(1,2))^2 );

% Compute pixel-to-meter conversion factor (pixels per meter)
pixelsPerMeter = rulerPixelLength / 0.45;

% Define ROI radius in pixels (0.04 m in pixels)
roiRadius_pixels = 0.04 * pixelsPerMeter;

% 2. Ask the user to select the centers for the two circular ROIs
title('Click on the center of the first ROI');
[center1_x, center1_y] = ginput(1);
center1 = [center1_x, center1_y];

title('Click on the center of the second ROI');
[center2_x, center2_y] = ginput(1);
center2 = [center2_x, center2_y];

% Set the same radius for both ROIs (in pixels)
radius1 = roiRadius_pixels;
radius2 = roiRadius_pixels;

% -------------------------------------------------------------------------
% Parameters for exploration detection
window_duration = 0.2; % Seconds (changed from 0.5)
min_overlap = 0.8; % Minimum fraction of frames within ROI to consider a valid window
movement_threshold = 1.5; % Pixels (lowered from 2)
% stable_frames_threshold will be set dynamically based on FPS *** Modification ***

% Initialize data structures for event logging
eventLogData = {}; % Cell array to store timestamps for each CSV

% Initialize summary data
summaryData = [];

% Loop through all CSV and corresponding video files
for i = 1:length(csv_files)
    % Load the CSV file
    csv_fileName = csv_files(i).name;
    fprintf('Reading CSV file: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);
    
    % Load the corresponding video file
    video_fileName = video_files(i).name;
    fprintf('Loading video file: %s\n', video_fileName);
    videoObj = VideoReader(video_fileName);
    
    % Retrieve FPS from the video file
    fps = videoObj.FrameRate;
    if isempty(fps) || isnan(fps)
        warning('FPS information missing in %s. Using default FPS = 30.', video_fileName);
        fps = 30; % Default FPS if not available
    end
    
    % *** Modification 1: Set stable_frames_threshold based on FPS ***
    stable_frames_threshold = 2 * fps; % 2 seconds based on actual FPS
    
    % *** Modification 2: Limit analysis to the first 6 minutes ***
    max_duration_sec = 6 * 60; % 6 minutes in seconds
    max_frame = round(max_duration_sec * fps); % Maximum frame number to analyze
    
    % Calculate window size based on window_duration and FPS
    window_size = round(window_duration * fps);
    if window_size < 1
        window_size = 1; % Ensure at least one frame
    end
    
    % Get the frame information from the CSV
    frames = csvData{:, 1};              % 1st column: Frame numbers
    noseX = csvData{:, 2};               % 2nd column: X coordinates of the nose
    noseY = csvData{:, 3};               % 3rd column: Y coordinates of the nose
    nose_likelihood = csvData{:, 4};     % 4th column: Likelihood of the nose detection
    headX = csvData{:, 5};               % 5th column: X coordinates of the head
    headY = csvData{:, 6};               % 6th column: Y coordinates of the head
    head_likelihood = csvData{:, 7};     % 7th column: Likelihood of the head detection
    bodyX = csvData{:, 8};               % 8th column: X coordinates of the body
    bodyY = csvData{:, 9};               % 9th column: Y coordinates of the body
    body_likelihood = csvData{:, 10};    % 10th column: Likelihood of the body detection
    
    % Exclude points with low likelihood
    validPoints = (nose_likelihood >= likelihood_threshold) & ...
                  (head_likelihood >= likelihood_threshold) & ...
                  (body_likelihood >= likelihood_threshold);
    noseX = noseX(validPoints);
    noseY = noseY(validPoints);
    headX = headX(validPoints);
    headY = headY(validPoints);
    bodyX = bodyX(validPoints);
    bodyY = bodyY(validPoints);
    frames_valid = frames(validPoints);
    
    % Exclude points outside the open field rectangle
    inOpenField = noseX >= openFieldPosition(1) & noseX <= (openFieldPosition(1) + openFieldPosition(3)) & ...
                  noseY >= openFieldPosition(2) & noseY <= (openFieldPosition(2) + openFieldPosition(4));
    noseX = noseX(inOpenField);
    noseY = noseY(inOpenField);
    headX = headX(inOpenField);
    headY = headY(inOpenField);
    bodyX = bodyX(inOpenField);
    bodyY = bodyY(inOpenField);
    frames_valid = frames_valid(inOpenField);
    
    % *** Modification 2: Limit to first 6 minutes ***
    % Find indices where frame number is less than or equal to max_frame
    within_time = frames_valid <= max_frame;
    noseX = noseX(within_time);
    noseY = noseY(within_time);
    headX = headX(within_time);
    headY = headY(within_time);
    bodyX = bodyX(within_time);
    bodyY = bodyY(within_time);
    frames_valid = frames_valid(within_time);
    
    % Calculate distances from nose to ROI centers (in pixels)
    dist1 = sqrt((noseX - center1(1)).^2 + (noseY - center1(2)).^2);
    dist2 = sqrt((noseX - center2(1)).^2 + (noseY - center2(2)).^2);
    
    % Determine if nose is within ROI (using the pixel-based radius)
    inROI1 = dist1 <= radius1;
    inROI2 = dist2 <= radius2;
    
    % Detect exploration for both ROIs
    [exploreROI1, exploreEventsROI1, totalTimeROI1, timestampsROI1] = detectExplorationBouts( ...
        inROI1, frames_valid, center1, radius1, bodyX, bodyY, noseX, noseY, fps, window_size, ...
        min_overlap, movement_threshold, stable_frames_threshold);

    [exploreROI2, exploreEventsROI2, totalTimeROI2, timestampsROI2] = detectExplorationBouts( ...
        inROI2, frames_valid, center2, radius2, bodyX, bodyY, noseX, noseY, fps, window_size, ...
        min_overlap, movement_threshold, stable_frames_threshold);


    % Log timestamps
    eventLogData{i,1} = csv_fileName;
    eventLogData{i,2} = timestampsROI1;
    eventLogData{i,3} = timestampsROI2;
    
    % Plot the trajectory on the first frame of this video
    videoObj.CurrentTime = 0; % Reset to first frame
    firstFrame = readFrame(videoObj);
    figure('Name', sprintf('Exploration Trajectory - %s', video_fileName), 'NumberTitle', 'off');
    imshow(firstFrame);
    hold on;
    
    % Plot all nose trajectories
    plot(noseX, noseY, 'm', 'LineWidth', 1); % Original nose trajectory
    
    % Highlight exploration trajectories
    exploreIndices1 = find(exploreROI1);
    exploreIndices2 = find(exploreROI2);
    plot(noseX(exploreIndices1), noseY(exploreIndices1), 'b.', 'MarkerSize', 10); % Exploration in ROI1
    plot(noseX(exploreIndices2), noseY(exploreIndices2), 'c.', 'MarkerSize', 10); % Exploration in ROI2
    
    % Draw the ROIs and open field
    viscircles(center1, radius1, 'EdgeColor', 'r'); % Draw the first circle
    viscircles(center2, radius2, 'EdgeColor', 'r'); % Draw the second circle
    rectangle('Position', openFieldPosition, 'EdgeColor', 'g', 'LineWidth', 2); % Draw the open field rectangle
    title(sprintf('Video: %s', video_fileName));
    legend('Nose Trajectory', 'Exploration ROI1', 'Exploration ROI2', 'ROIs', 'Location', 'Best');
    hold off;
    
    % Save the frame with circles, trajectory, and exploration as a TIFF file
    [~, videoFileNameNoExt, ~] = fileparts(video_fileName); % Remove extension from video file name
    saveas(gcf, sprintf('%s_FirstFrame_with_ROIs_and_Exploration.tif', videoFileNameNoExt));
    
    % Append the results to the summary data
    summaryData = [summaryData; {csv_fileName, totalTimeROI1, exploreEventsROI1, totalTimeROI2, exploreEventsROI2}];
    
    % Display the results for this video
    fprintf('Video %d - Time spent exploring object 1: %.2f seconds with %d events\n', i, totalTimeROI1, exploreEventsROI1);
    fprintf('Video %d - Time spent exploring object 2: %.2f seconds with %d events\n', i, totalTimeROI2, exploreEventsROI2);
end

% Save the summary data to an Excel file (with discrimination ratios + notes)
summaryTable = cell2table(summaryData, 'VariableNames', ...
    {'CSV_File', 'Time_Object_1_sec', 'Events_Object_1', 'Time_Object_2_sec', 'Events_Object_2'});

t1 = summaryTable.Time_Object_1_sec;
t2 = summaryTable.Time_Object_2_sec;
den = t1 + t2;

% Discrimination ratios
DR_1vs2 = nan(height(summaryTable),1);
DR_2vs1 = nan(height(summaryTable),1);

validDen = den > 0; % avoid divide-by-zero
DR_1vs2(validDen) = (t1(validDen) - t2(validDen)) ./ den(validDen);
DR_2vs1(validDen) = (t2(validDen) - t1(validDen)) ./ den(validDen);

summaryTable.DR_Object1_vs_Object2 = DR_1vs2;
summaryTable.DR_Object2_vs_Object1 = DR_2vs1;

% Notes / warnings
notes = strings(height(summaryTable),1);

% Rule 1: any zero in either object time
zeroMask = (t1 == 0) | (t2 == 0);
notes(zeroMask) = notes(zeroMask) + "WARNING: zero exploring time detected (Time_Object_1_sec or Time_Object_2_sec = 0). ";

% Rule 2: total exploring time < 5 sec
lowTotalMask = den < 5;
notes(lowTotalMask) = notes(lowTotalMask) + "WARNING: total exploration time < 5 s (possible low movement). Check tracking trace. ";

% If denominator = 0, ratios are undefined
denZeroMask = den == 0;
notes(denZeroMask) = notes(denZeroMask) + "WARNING: Time_Object_1_sec + Time_Object_2_sec = 0, discrimination ratio undefined. ";

summaryTable.Note = strtrim(notes);

writetable(summaryTable, 'summary_NPR.xlsx', 'Sheet', 'Summary');

% Prepare event log table with aligned timestamps
eventLogStruct = struct();
for i = 1:length(csv_files)
    eventLogStruct.CSV_File{i,1} = eventLogData{i,1};
    
    % Format ROI1 timestamps as comma-separated strings
    if ~isempty(eventLogData{i,2})
        str_ROI1 = sprintf('%.2f, ', eventLogData{i,2});
        str_ROI1 = strtrim(str_ROI1);
        str_ROI1(end) = []; % Remove trailing comma
    else
        str_ROI1 = '';
    end
    eventLogStruct.Exploration_ROI1_Timestamps{i,1} = str_ROI1;
    
    % Format ROI2 timestamps as comma-separated strings
    if ~isempty(eventLogData{i,3})
        str_ROI2 = sprintf('%.2f, ', eventLogData{i,3});
        str_ROI2 = strtrim(str_ROI2);
        str_ROI2(end) = []; % Remove trailing comma
    else
        str_ROI2 = '';
    end
    eventLogStruct.Exploration_ROI2_Timestamps{i,1} = str_ROI2;
end
eventLogTable = struct2table(eventLogStruct);
writetable(eventLogTable, 'summary_NPR.xlsx', 'Sheet', 'Event_Logs');

fprintf('Summary data saved to summary_NPR.xlsx\n');

% ------------------------- Helper Functions -------------------------

function [exploreFlags, numEvents, totalTime, timestamps, boutTable] = detectExplorationBouts( ...
    inROI, frames_valid, center, radius, bodyX, bodyY, noseX, noseY, fps, window_size, min_overlap, ...
    movement_threshold, stable_frames_threshold)

    % ---- Tunable bout logic ----
    merge_gap_sec  = 0.40;  % merge bouts if gap <= this (sec). Try 0.2–0.6
    min_bout_sec   = 0.15;  % discard bouts shorter than this (sec). Try 0.15–0.5
    % ----------------------------

    N = numel(inROI);
    exploreFlags = false(N,1);
    timestamps   = [];
    boutTable    = table();

    if N == 0
        numEvents = 0; totalTime = 0;
        return;
    end

    % 1) Windowed ROI-occupancy criterion -> candidate frames
    min_frames_in_ROI = ceil(min_overlap * window_size);

    % movsum gives, for each frame, how many of the last window_size frames were in ROI
    roiCount = movsum(double(inROI), [window_size-1, 0]);  % trailing window
    candidate = roiCount >= min_frames_in_ROI;

    % 2) Optional: require "approach" (distance trending down) to reduce false positives
    dist = hypot(noseX - center(1), noseY - center(2));
    dDist = [0; diff(dist)];
    approach = movmean(dDist, [window_size-1, 0]) < 0; % negative trend
    candidate = candidate & approach;

    % 3) Exclude "sitting on object" using vectorized stability test
    bodyDist = hypot(bodyX - center(1), bodyY - center(2));
    overlap = bodyDist <= (radius/2);

    % stability: small movement of body over stable_frames_threshold
    stableX = movstd(bodyX, [stable_frames_threshold-1, 0], 'omitnan') < movement_threshold;
    stableY = movstd(bodyY, [stable_frames_threshold-1, 0], 'omitnan') < movement_threshold;
    stable = stableX & stableY;

    sitting = overlap & stable;

    exploreRaw = candidate & ~sitting;

    % 4) Merge bouts separated by brief gaps (tracking jitter / brief exits)
    merge_gap_frames = max(0, round(merge_gap_sec * fps));
    exploreMerged = mergeLogicalBouts(exploreRaw, merge_gap_frames);

    % 5) Remove very short bouts
    min_bout_frames = max(1, round(min_bout_sec * fps));
    exploreFlags = removeShortBouts(exploreMerged, min_bout_frames);

    % 6) Extract bout start/end indices
    [starts, ends] = logicalBouts(exploreFlags);

    numEvents = numel(starts);
    if numEvents == 0
        totalTime = 0;
        timestamps = [];
        boutTable = table();
        return;
    end

    durations_sec = (ends - starts + 1) / fps;
    totalTime = sum(durations_sec);

    % Use true timestamps based on frames_valid (better than j/fps after filtering)
    timestamps = frames_valid(starts) / fps;

    boutTable = table( ...
        frames_valid(starts), frames_valid(ends), timestamps, durations_sec, ...
        'VariableNames', {'StartFrame','EndFrame','StartTimeSec','DurationSec'});
end

function out = mergeLogicalBouts(x, maxGap)
    % Merge bouts separated by <= maxGap false frames
    out = x(:);
    if maxGap <= 0
        return;
    end
    [s,e] = logicalBouts(out);
    if isempty(s), return; end

    for k = 1:numel(s)-1
        gap = s(k+1) - e(k) - 1;
        if gap > 0 && gap <= maxGap
            out(e(k)+1 : s(k+1)-1) = true;
        end
    end
end

function out = removeShortBouts(x, minLen)
    out = x(:);
    [s,e] = logicalBouts(out);
    for k = 1:numel(s)
        if (e(k) - s(k) + 1) < minLen
            out(s(k):e(k)) = false;
        end
    end
end

function [starts, ends] = logicalBouts(x)
    x = x(:);
    d = diff([false; x; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
end
