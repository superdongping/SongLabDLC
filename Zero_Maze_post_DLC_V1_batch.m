clear all;
close all;
clc;

% Set the likelihood threshold for the open arms
likelihood_threshold = 0.3;

% Search for CSV and MP4 files in the current folder
csv_files = dir('*.csv');
video_files = dir('*.mp4');

if isempty(csv_files) || isempty(video_files)
    error('No CSV or MP4 files found in the current directory.');
end

% Define constants
frames_per_second = 30;  % Video frame rate (30 FPS)
analysis_duration_sec = 6 * 60;  % Analyze first 6 minutes (360 seconds)
max_frames_to_analyze = frames_per_second * analysis_duration_sec;  % Maximum number of frames to analyze

% First pass: find the global maximum for consistent color scaling
max_time_spent = 0;

% Initialize summary data
summaryData = [];

% Loop through all videos and CSV files to find the maximum time spent in any bin
for i = 1:length(csv_files)
    % Load the CSV file
    csv_fileName = csv_files(i).name;
    fprintf('Reading CSV file: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);
    
    % Load the video file to get dimensions
    video_fileName = video_files(i).name;
    fprintf('Loading video file: %s\n', video_fileName);
    videoObj = VideoReader(video_fileName);
    
    % Get the frame information from the CSV
    frames = csvData{:, 1};      % 1st column: Frame numbers
    MainBodyX = csvData{:, 8};   % X coordinates of the MainBody
    MainBodyY = csvData{:, 9};   % Y coordinates of the MainBody
    likelihood = csvData{:, 10}; % Likelihood/confidence of the detection
    
    % Limit the analysis to the first 6 minutes (or available frames)
    num_frames = min(max_frames_to_analyze, length(frames));
    frames = frames(1:num_frames);
    MainBodyX = MainBodyX(1:num_frames);
    MainBodyY = MainBodyY(1:num_frames);
    likelihood = likelihood(1:num_frames);
    
    % Determine which points are in the open arms (likelihood >= threshold)
    validPoints = likelihood >= likelihood_threshold;
    
    % Create a grid for the heatmap using video frame dimensions
    x_edges = linspace(min(MainBodyX), max(MainBodyX), 100); % 100 bins along X-axis
    y_edges = linspace(min(MainBodyY), max(MainBodyY), 100); % 100 bins along Y-axis

    % Accumulate the time spent in each bin (i.e., coordinate)
    timeSpent = histcounts2(MainBodyX(validPoints), MainBodyY(validPoints), x_edges, y_edges);
    
    % Track the maximum time spent in any bin across all videos
    max_time_spent = max(max_time_spent, max(timeSpent(:)));
end

% Second pass: plot the heatmaps using the global max time for color scaling
for i = 1:length(csv_files)
    % Load the CSV file
    csv_fileName = csv_files(i).name;
    fprintf('Reading CSV file: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);
    
    % Load the video file to get dimensions
    video_fileName = video_files(i).name;
    fprintf('Loading video file: %s\n', video_fileName);
    videoObj = VideoReader(video_fileName);
    firstFrame = readFrame(videoObj); % Get the first frame for plotting
    
    % Get the frame information from the CSV
    frames = csvData{:, 1};      % 1st column: Frame numbers
    MainBodyX = csvData{:, 8};   % X coordinates of the MainBody
    MainBodyY = csvData{:, 9};   % Y coordinates of the MainBody
    likelihood = csvData{:, 10}; % Likelihood/confidence of the detection
    
    % Limit the analysis to the first 6 minutes (or available frames)
    num_frames = min(max_frames_to_analyze, length(frames));
    frames = frames(1:num_frames);
    MainBodyX = MainBodyX(1:num_frames);
    MainBodyY = MainBodyY(1:num_frames);
    likelihood = likelihood(1:num_frames);
    
    % Determine which points are in the open arms (likelihood >= threshold)
    validPoints = likelihood >= likelihood_threshold;

    % Create a grid for the heatmap using video frame dimensions
    x_edges = linspace(min(MainBodyX), max(MainBodyX), 100); % 100 bins along X-axis
    y_edges = linspace(min(MainBodyY), max(MainBodyY), 100); % 100 bins along Y-axis

    % Accumulate the time spent in each bin (i.e., coordinate)
    timeSpent = histcounts2(MainBodyX(validPoints), MainBodyY(validPoints), x_edges, y_edges);
    
    % Normalize the time spent to the global max
    upper_limit = 0.8 * max_time_spent; % Adjust the upper limit (80% of max for better contrast)
    timeSpentNormalized = timeSpent / upper_limit; 
    
    % ---- Apply Gaussian Smoothing ---- %
    sigma = 2;  % Set the standard deviation for Gaussian smoothing
    timeSpentSmoothed = imgaussfilt(timeSpentNormalized, sigma);  % Apply Gaussian filter
    
    % ---- Plot the valid MainBody points on the first frame ---- %
    figure;
    imshow(firstFrame); % Display the first frame
    hold on;
    scatter(MainBodyX(validPoints), MainBodyY(validPoints), 8, 'g', 'filled', 'MarkerFaceAlpha', 0.5);  % Plot valid points
    title(sprintf('Video: %s - Open Arm Exploration', video_fileName));
    hold off;
    
    % Save the frame with the points as a TIFF file
    [~, videoFileNameNoExt, ~] = fileparts(video_fileName);
    saveas(gcf, sprintf('%s_OpenArms_FirstFrame.tif', videoFileNameNoExt));

    % ---- Smoothed Heatmap Generation (without frame overlay) ---- %
    figure;
    imagesc(x_edges, y_edges, timeSpentSmoothed'); % Plot the smoothed heatmap (transposed for correct orientation)
    colormap(jet); % Use the 'jet' colormap for blue-to-red gradient
    colorbar; % Add a colorbar to show the scale
    caxis([0 0.05]); % Set colorbar scale consistent across all mice (based on the normalized scale)
    title('Mouse Heatmap - Time Spent in Each Area (Smoothed)');
    
    % Make the size of the heatmap consistent with the video frame
    axis equal;
    set(gca, 'XLim', [min(MainBodyX) max(MainBodyX)], 'YLim', [min(MainBodyY) max(MainBodyY)]);
    
    % Save the smoothed heatmap as a TIFF file
    saveas(gcf, sprintf('%s_Heatmap_jet_Smoothed.tif', videoFileNameNoExt));
    
    % Calculate the total trajectory distance
    total_distance = sum(sqrt(diff(MainBodyX(validPoints)).^2 + diff(MainBodyY(validPoints)).^2));  % Total distance in pixels
    
    % Calculate time spent in open arms and ratio
    time_in_open_arms = sum(validPoints) / frames_per_second;  % Time spent in open arms in seconds
    total_time = num_frames / frames_per_second;  % Total analysis time in seconds
    ratio_open_arms = time_in_open_arms / total_time;  % Ratio of time spent in open arms
    
    % Append the results to the summary data
    summaryData = [summaryData; {csv_fileName, total_distance, time_in_open_arms, ratio_open_arms}];
    
    % Display results for this video
    fprintf('Video %d - Total distance: %.2f pixels, Time in open arms: %.2f sec, Ratio: %.2f\n', ...
            i, total_distance, time_in_open_arms, ratio_open_arms);
end

% Save the summary data to an Excel file
summaryTable = cell2table(summaryData, 'VariableNames', {'CSV_File', 'Total_Trajectory_Distance_pixels', 'Time_in_Open_Arms_sec', 'Ratio_Open_Arms'});
writetable(summaryTable, 'summary_Zero_maze.xlsx');
fprintf('Summary data saved to summary_Zero_maze.xlsx\n');
