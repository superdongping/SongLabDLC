clear all; 
close all;
warning off;
clc;

% Set the likelihood threshold
likelihood_threshold = 0.05;

% Search for CSV and MP4 files in the current folder
csv_files = dir('*.csv');
video_files = dir('*.mp4');

if isempty(csv_files) || isempty(video_files)
    error('No CSV or MP4 files found in the current directory.');
end

% Load the first video to define the analysis area
video_fileName = video_files(1).name; 
fprintf('Loading first video: %s\n', video_fileName);

% Load the first video and get the first frame
videoObj = VideoReader(video_fileName);
firstFrame = readFrame(videoObj);

% Display the first frame and let the user draw the outer region for analysis
figure;
imshow(firstFrame);
title('Draw a rectangle to define the overall region for analysis');
outerRegionRect = drawrectangle('Color', 'c'); % User defines the outer region
outerRegionPosition = outerRegionRect.Position; % [x, y, width, height]

% Draw the open field box (bottom layer)
title('Draw a rectangle to define the open field area');
openFieldRect = drawrectangle('Color', 'g'); % User defines the open field area
openFieldPosition = openFieldRect.Position; % [x, y, width, height]

% Let the user draw a line to define 0.45m scale on the open field
title('Draw a line to define 0.45m scale on one side of the open field');
scaleLine = drawline('Color','b'); 
scaleLengthInPixels = sqrt(sum((scaleLine.Position(2, :) - scaleLine.Position(1, :)).^2));
scaleInMeters = 0.45; % Actual length in meters
pixelsPerMeter = scaleLengthInPixels / scaleInMeters; % Conversion factor

% Define the central area (0.25 m square in the center of the open field)
centralWidthPixels = 0.25 * pixelsPerMeter;
centralAreaX = openFieldPosition(1) + (openFieldPosition(3) - centralWidthPixels) / 2;
centralAreaY = openFieldPosition(2) + (openFieldPosition(4) - centralWidthPixels) / 2;
centralArea = [centralAreaX, centralAreaY, centralWidthPixels, centralWidthPixels];

% Initialize variables for heatmap normalization
summaryData = [];
heatmapMin = inf;
heatmapMax = -inf;
heatmapDimensions = [100, 100]; % Define a fixed size for the heatmap grid

% First pass: calculate min and max values for the heatmap normalization
for i = 1:length(csv_files)
    % Load the CSV file
    csv_fileName = csv_files(i).name;
    fprintf('Reading CSV file: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);
    
    % Load the video file
    video_fileName = video_files(i).name;
    fprintf('Loading video file: %s\n', video_fileName);
    videoObj = VideoReader(video_fileName);
    [~, videoFileNameNoExt, ~] = fileparts(video_fileName); % Extract the file name without extension

    % Get the frame information from the CSV
    MainBodyX = csvData{:, 8};   % X coordinates of the MainBody
    MainBodyY = csvData{:, 9};   % Y coordinates of the MainBody
    likelihood = csvData{:, 10}; % Likelihood/confidence of the detection
    
    % Exclude points with low likelihood
    validPoints = likelihood >= likelihood_threshold;
    MainBodyX = MainBodyX(validPoints);
    MainBodyY = MainBodyY(validPoints);
    
    % Exclude points outside the outer region
    inOuterRegion = MainBodyX >= outerRegionPosition(1) & MainBodyX <= (outerRegionPosition(1) + outerRegionPosition(3)) & ...
                    MainBodyY >= outerRegionPosition(2) & MainBodyY <= (outerRegionPosition(2) + outerRegionPosition(4));
    MainBodyX = MainBodyX(inOuterRegion);
    MainBodyY = MainBodyY(inOuterRegion);
    
    % Exclude points outside the open field rectangle
    inOpenField = MainBodyX >= openFieldPosition(1) & MainBodyX <= (openFieldPosition(1) + openFieldPosition(3)) & ...
                  MainBodyY >= openFieldPosition(2) & MainBodyY <= (openFieldPosition(2) + openFieldPosition(4));
    MainBodyX = MainBodyX(inOpenField);
    MainBodyY = MainBodyY(inOpenField);
    
    % Generate heatmap of the MainBody positions
    heatmapData = histcounts2(MainBodyX, MainBodyY, ...
        linspace(openFieldPosition(1), openFieldPosition(1) + openFieldPosition(3), heatmapDimensions(1)), ...
        linspace(openFieldPosition(2), openFieldPosition(2) + openFieldPosition(4), heatmapDimensions(2)), ...
        'Normalization', 'probability');
    
    % Apply Gaussian smoothing to the heatmap
    sigma = 2;  % Set the standard deviation for Gaussian smoothing
    timeSpentSmoothed = imgaussfilt(heatmapData, sigma);

    % Update global min and max for consistent scaling
    heatmapMin = min(heatmapMin, min(timeSpentSmoothed(:)));
    heatmapMax = max(heatmapMax, max(timeSpentSmoothed(:)));
end

% Adjust the maximum value for brightness enhancement
adjustedHeatmapMax = heatmapMax * 0.20;  % Set the heatmap display max to 10% of the computed max

% Second pass: process each video and save results with consistent heatmap scaling
for i = 1:length(csv_files)
    % Load the CSV file
    csv_fileName = csv_files(i).name;
    fprintf('Reading CSV file: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);
    
    % Load the video file
    video_fileName = video_files(i).name;
    fprintf('Loading video file: %s\n', video_fileName);
    videoObj = VideoReader(video_fileName);
    [~, videoFileNameNoExt, ~] = fileparts(video_fileName);

    % Get the frame information from the CSV
    MainBodyX = csvData{:, 8};
    MainBodyY = csvData{:, 9};
    likelihood = csvData{:, 10};
    
    % Exclude points with low likelihood
    validPoints = likelihood >= likelihood_threshold;
    MainBodyX = MainBodyX(validPoints);
    MainBodyY = MainBodyY(validPoints);
    
    % Exclude points outside the outer region
    inOuterRegion = MainBodyX >= outerRegionPosition(1) & MainBodyX <= (outerRegionPosition(1) + outerRegionPosition(3)) & ...
                    MainBodyY >= outerRegionPosition(2) & MainBodyY <= (outerRegionPosition(2) + outerRegionPosition(4));
    MainBodyX = MainBodyX(inOuterRegion);
    MainBodyY = MainBodyY(inOuterRegion);
    
    % Exclude points outside the open field rectangle
    inOpenField = MainBodyX >= openFieldPosition(1) & MainBodyX <= (openFieldPosition(1) + openFieldPosition(3)) & ...
                  MainBodyY >= openFieldPosition(2) & MainBodyY <= (openFieldPosition(2) + openFieldPosition(4));
    MainBodyX = MainBodyX(inOpenField);
    MainBodyY = MainBodyY(inOpenField);
    
    % Calculate the total trajectory distance
    trajectoryDistance = sum(sqrt(diff(MainBodyX).^2 + diff(MainBodyY).^2)) / pixelsPerMeter;
    
    % Calculate time spent in the central area
    inCentralArea = MainBodyX >= centralArea(1) & MainBodyX <= (centralArea(1) + centralArea(3)) & ...
                    MainBodyY >= centralArea(2) & MainBodyY <= (centralArea(2) + centralArea(4));
    timeInCentralArea = sum(inCentralArea) / 30;
    totalTime = length(MainBodyX) / 30;
    centralTimeRatio = timeInCentralArea / totalTime;
    
    % Calculate distance in the central area and speeds
    centralX = MainBodyX(inCentralArea);
    centralY = MainBodyY(inCentralArea);
    distanceInCentral = sum(sqrt(diff(centralX).^2 + diff(centralY).^2)) / pixelsPerMeter;
    speedTotal = trajectoryDistance / totalTime;
    speedCentral = distanceInCentral / timeInCentralArea;
    
    % Generate and smooth heatmap with consistent size
    heatmapData = histcounts2(MainBodyX, MainBodyY, ...
        linspace(openFieldPosition(1), openFieldPosition(1) + openFieldPosition(3), heatmapDimensions(1)), ...
        linspace(openFieldPosition(2), openFieldPosition(2) + openFieldPosition(4), heatmapDimensions(2)), ...
        'Normalization', 'probability');
    timeSpentSmoothed = imgaussfilt(heatmapData, sigma);
    
    % Plot and save the smoothed heatmap with consistent scaling
    figure;
    
    % Flip the heatmap data along the Y-axis to match the video orientation
    imagesc(flipud(timeSpentSmoothed), [heatmapMin, adjustedHeatmapMax]);
    colormap('jet');
    colorbar;
    title('Smoothed Heatmap of MainBody Trajectory');
    xlabel('X Position');
    ylabel('Y Position');
    set(gca, 'YDir', 'normal');
    
    saveas(gcf, sprintf('%s_SmoothedHeatmap_Adjusted.tif', videoFileNameNoExt));
    close(gcf);

    % Save trajectory plot
    firstFrame = readFrame(videoObj);
    figure;
    imshow(firstFrame);
    hold on;
    plot(MainBodyX, MainBodyY, 'm', 'LineWidth', 1);
    rectangle('Position', openFieldPosition, 'EdgeColor', 'g', 'LineWidth', 2);
    rectangle('Position', centralArea, 'EdgeColor', 'r', 'LineWidth', 2);
    line(scaleLine.Position(:, 1), scaleLine.Position(:, 2), 'Color', 'b', 'LineWidth', 2);
    title(sprintf('Video: %s', video_fileName));
    hold off;
    saveas(gcf, sprintf('%s_FirstFrame_with_Trajectory.tif', videoFileNameNoExt));
    close(gcf);
    
    % Append results to summary data
    summaryData = [summaryData; {csv_fileName, trajectoryDistance, timeInCentralArea, centralTimeRatio, ...
        distanceInCentral, speedTotal, speedCentral}];
end

% Save the summary data to an Excel file
summaryTable = cell2table(summaryData, 'VariableNames', {'CSV_File', 'Total_Trajectory_Distance_m', ...
    'Time_in_Central_s', 'Central_Time_Ratio', 'Distance_in_Central_m', 'Speed_Total_mps', 'Speed_Central_mps'});
writetable(summaryTable, 'summary_Open_field.xlsx');
fprintf('Summary data saved to summary_Open_field.xlsx\n');
