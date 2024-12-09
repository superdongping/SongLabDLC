clear all;
close all;
clc;

% Set the likelihood threshold
% (If you still need to filter points by likelihood, keep this. Otherwise, set it very low.)
likelihood_threshold = 0.01;

% Search for CSV and MP4 files in the current folder
csv_files = dir('*.csv');
video_files = dir('*.mp4');

if isempty(csv_files) || isempty(video_files)
    error('No CSV or MP4 files found in the current directory.');
end

% --- Obtain FPS from the first video and round it up ---
firstVideoObj = VideoReader(video_files(1).name);
frames_per_second = ceil(firstVideoObj.FrameRate);

% Analyze first 6 minutes
analysis_duration_sec = 6 * 60;  % 360 seconds
max_frames_to_analyze = frames_per_second * analysis_duration_sec;

% --- Ask user to define closed arm ROIs using the first frame of the first video ---
firstFrame = readFrame(firstVideoObj);
figure('Name','Draw Closed Arm ROIs');
imshow(firstFrame);
title('Draw the polygon for the first closed arm ROI, double-click to finish');
closedArmsROI1 = drawpolygon;
wait(closedArmsROI1); % Wait until user finishes drawing

title('Draw the polygon for the second closed arm ROI, double-click to finish');
closedArmsROI2 = drawpolygon;
wait(closedArmsROI2);

maskROI1 = closedArmsROI1.Position;
maskROI2 = closedArmsROI2.Position;
close(gcf);

% First pass: find global maximum for heatmap scaling
max_time_spent = 0;

% Initialize summary data
summaryData = [];

% Loop through all videos/CSVs to find global max
for i = 1:length(csv_files)
    csv_fileName = csv_files(i).name;
    video_fileName = video_files(i).name;
    fprintf('Reading CSV file: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);

    fprintf('Loading video file: %s\n', video_fileName);
    videoObj = VideoReader(video_fileName);

    % Extract tracked coordinates
    frames = csvData{:, 1};      % Frame numbers
    MainBodyX = csvData{:, 8};   % X coordinates of the MainBody
    MainBodyY = csvData{:, 9};   % Y coordinates of the MainBody
    likelihood = csvData{:, 10}; % Likelihood/confidence of the detection

    % Limit the analysis to first 6 mins (or available frames)
    num_frames = min(max_frames_to_analyze, length(frames));
    frames = frames(1:num_frames);
    MainBodyX = MainBodyX(1:num_frames);
    MainBodyY = MainBodyY(1:num_frames);
    likelihood = likelihood(1:num_frames);

    % Valid points based on likelihood threshold
    validPoints = likelihood >= likelihood_threshold;

    % Compute edges for heatmap
    x_edges = linspace(min(MainBodyX), max(MainBodyX), 100);
    y_edges = linspace(min(MainBodyY), max(MainBodyY), 100);

    % Time spent in each bin
    timeSpent = histcounts2(MainBodyX(validPoints), MainBodyY(validPoints), x_edges, y_edges);
    max_time_spent = max(max_time_spent, max(timeSpent(:)));
end

% Second pass: generate plots and compute metrics
for i = 1:length(csv_files)
    csv_fileName = csv_files(i).name;
    video_fileName = video_files(i).name;
    fprintf('Processing: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);

    videoObj = VideoReader(video_fileName);
    firstFrame = readFrame(videoObj); % First frame of this video

    frames = csvData{:, 1};
    MainBodyX = csvData{:, 8};
    MainBodyY = csvData{:, 9};
    likelihood = csvData{:, 10};

    % Limit to first 6 minutes
    num_frames = min(max_frames_to_analyze, length(frames));
    frames = frames(1:num_frames);
    MainBodyX = MainBodyX(1:num_frames);
    MainBodyY = MainBodyY(1:num_frames);
    likelihood = likelihood(1:num_frames);

    validPoints = likelihood >= likelihood_threshold;

    % Create heatmap edges
    x_edges = linspace(min(MainBodyX), max(MainBodyX), 100);
    y_edges = linspace(min(MainBodyY), max(MainBodyY), 100);

    % Time spent
    timeSpent = histcounts2(MainBodyX(validPoints), MainBodyY(validPoints), x_edges, y_edges);

    % Normalize timeSpent
    upper_limit = 0.8 * max_time_spent;
    timeSpentNormalized = timeSpent / upper_limit;

    % Apply Gaussian smoothing
    sigma = 2;
    timeSpentSmoothed = imgaussfilt(timeSpentNormalized, sigma);

    % ----- Determine open vs closed arm time -----
    % Convert each point coordinate into a decision if it's in closed arm or not
    % For each frame, if valid, check if (X,Y) is inside either closed arm ROI
    inClosedArm1 = inpolygon(MainBodyX, MainBodyY, maskROI1(:,1), maskROI1(:,2));
    inClosedArm2 = inpolygon(MainBodyX, MainBodyY, maskROI2(:,1), maskROI2(:,2));

    inClosedArm = (inClosedArm1 | inClosedArm2) & validPoints; % Valid points inside closed arms
    closed_arm_frames = sum(inClosedArm);  % number of frames inside closed arms
    total_time = num_frames / frames_per_second;
    time_in_closed_arms = closed_arm_frames / frames_per_second;

    % The rest is open arms time
    time_in_open_arms = total_time - time_in_closed_arms;
    ratio_open_arms = time_in_open_arms / total_time;

    % Calculate total distance traveled (only between valid consecutive points)
    validX = MainBodyX(validPoints);
    validY = MainBodyY(validPoints);
    total_distance = sum(sqrt(diff(validX).^2 + diff(validY).^2));

    % ---- Plot trajectory on first frame ----
    figure;
    imshow(firstFrame);
    hold on;
    scatter(MainBodyX(validPoints), MainBodyY(validPoints), 8, 'g', 'filled', 'MarkerFaceAlpha', 0.5);
    title(sprintf('Video: %s - Trajectory', video_fileName), 'Interpreter', 'none');
    hold off;

    % Save the frame with trajectory as TIFF
    [~, videoFileNameNoExt, ~] = fileparts(video_fileName);
    saveas(gcf, sprintf('%s_Trajectory_FirstFrame.tif', videoFileNameNoExt));

    % ---- Plot Smoothed Heatmap ----
    figure;
    imagesc(x_edges, y_edges, timeSpentSmoothed');
    colormap(jet);
    colorbar;
    caxis([0 0.05]);
    title(sprintf('Mouse Heatmap - %s', video_fileName), 'Interpreter', 'none');
    axis equal;
    set(gca, 'XLim', [min(MainBodyX) max(MainBodyX)], 'YLim', [min(MainBodyY) max(MainBodyY)]);

    % Save smoothed heatmap
    saveas(gcf, sprintf('%s_Heatmap_jet_Smoothed.tif', videoFileNameNoExt));

    % Append summary data
    summaryData = [summaryData; {csv_fileName, total_distance, time_in_open_arms, ratio_open_arms}];

    fprintf('Video %d - Total distance: %.2f pixels, Time in open arms: %.2f sec, Ratio open arms: %.2f\n', ...
        i, total_distance, time_in_open_arms, ratio_open_arms);
end

% Save the summary to Excel
summaryTable = cell2table(summaryData, 'VariableNames', {'CSV_File','Total_Distance_pixels','Time_in_Open_Arms_sec','Ratio_Open_Arms'});
writetable(summaryTable, 'summary_Zero_maze.xlsx');
fprintf('Summary data saved to summary_Zero_maze.xlsx\n');
