% MATLAB Script to Quantify and Visualize Weighted Immobility in Tail Suspension Test for Multiple Mice
% This script processes multiple CSV files in the current folder, uses the first MP4 file for scaling,
% retrieves FPS from the first MP4 file, performs analysis from 0 to 360 seconds, saves the results accordingly,
% and generates annotated output videos with immobility status indicators.

% Clear workspace, close all figures, and clear command window
close all;
clear; clc;

%% Step 1: Get All CSV and MP4 Files in Current Directory

% Get list of all CSV files
csvFiles = dir('*.csv');
if isempty(csvFiles)
    error('No CSV files found in the current directory.');
end

% Get list of all MP4 files
mp4Files = dir('*.mp4');
if isempty(mp4Files)
    error('No MP4 video files found in the current directory. At least one MP4 file is required for scaling.');
end

% Use the first MP4 file for scaling and FPS
scalingMP4 = mp4Files(1);
videoFileFullPath = fullfile(scalingMP4.folder, scalingMP4.name);
fprintf('Using MP4 file for scaling and FPS: %s\n', scalingMP4.name);

%% Step 2: Load the First Frame and Retrieve FPS from the First MP4 File

% Create a VideoReader object for the scaling MP4
vidReaderScaling = VideoReader(videoFileFullPath);

% Retrieve FPS from the scaling MP4 file and round it to the nearest integer
fps_original = vidReaderScaling.FrameRate;
fps = round(fps_original);
fprintf('Original FPS of scaling video: %.2f\n', fps_original);
fprintf('Rounded FPS to be used: %d\n', fps);

% Read the first frame from the video
firstFrame = readFrame(vidReaderScaling);

% Display the first frame and let the user draw a line to indicate beam width
figure('Name', 'Select Beam Width for Scaling', 'NumberTitle', 'off');
imshow(firstFrame);
title('Draw a line to indicate the width of the beam (52 cm)');
h = drawline('Color', 'red', 'LineWidth', 2);
wait(h); % Wait until the user finishes drawing

% Get the position of the drawn line
linePos = h.Position; % [x1 y1; x2 y2]
pixelLength = sqrt((linePos(2,1) - linePos(1,1))^2 + (linePos(2,2) - linePos(1,2))^2);
knownWidth_cm = 52; % Known beam width in cm
pixels_to_cm = knownWidth_cm / pixelLength;
fprintf('Pixel to cm scaling factor: %.4f cm/pixel\n', pixels_to_cm);

% Close the figure after selection
close;

%% Step 3: Initialize Summary Data

% Initialize summary data as a cell array
summaryData = {};

%% Step 4: Loop Through Each CSV File (Each Mouse)

for i = 1:length(csvFiles)
    fprintf('\nProcessing Mouse %d/%d: %s\n', i, length(csvFiles), csvFiles(i).name);
    
    %% Step 4.1: Read and Parse the CSV Data
    
    csvFileFullPath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    fid = fopen(csvFileFullPath, 'r');
    if fid == -1
        warning('Cannot open the CSV file: %s. Skipping this mouse.', csvFileFullPath);
        continue;
    end
    
    % Read the first three header lines
    scorerLine = fgetl(fid);
    bodyPartsLine = fgetl(fid);
    coordsLine = fgetl(fid);
    
    % Split the header lines by comma
    scorerHeaders = strsplit(scorerLine, ',');
    bodyPartsHeaders = strsplit(bodyPartsLine, ',');
    coordsHeaders = strsplit(coordsLine, ',');
    
    % Verify that the number of headers align
    numColumns = length(scorerHeaders);
    if length(bodyPartsHeaders) ~= numColumns || length(coordsHeaders) ~= numColumns
        fclose(fid);
        warning('Mismatch in the number of columns among the header lines for %s. Skipping.', csvFiles(i).name);
        continue;
    end
    
    % Extract body part names and coordinate types (excluding the first column which is 'scorer')
    bodyParts = bodyPartsHeaders(2:end); % Cell array (e.g., {'Head', 'Head', 'Head', 'MainBody', ...})
    coordTypes = coordsHeaders(2:end);    % Should be {'x', 'y', 'likelihood', 'x', 'y', 'likelihood', ...}
    
    % Determine the number of body parts
    % Each body part has three columns: x, y, likelihood
    if mod((numColumns -1), 3) ~= 0
        fclose(fid);
        warning('Each body part should have exactly three columns: x, y, likelihood for %s. Skipping.', csvFiles(i).name);
        continue;
    end
    numBodyParts = (numColumns -1) / 3;
    
    % Assign unique identifiers to each body part occurrence
    % This is important if body parts are repeated
    bodyPartIdentifiers = cell(1, numBodyParts);
    for bp = 1:numBodyParts
        % Assuming bodyParts are grouped as x, y, likelihood
        bodyPartName = bodyParts{(bp-1)*3 + 1}; 
        bodyPartIdentifiers{bp} = sprintf('%s_%d', bodyPartName, bp);
    end
    
    % Initialize a cell array to store coordinates for each frame
    bodyCoordinates = {}; % To be filled with structures for each frame
    
    % Read the data starting from the fourth line
    % Initialize frame index
    frameIdx = 0;
    
    while ~feof(fid)
        dataLine = fgetl(fid);
        if ischar(dataLine)
            frameIdx = frameIdx + 1;
            data = strsplit(dataLine, ',');
            
            % Convert data to numeric, handling possible non-numeric entries
            numericData = str2double(data);
            
            if isnan(numericData(1))
                warning('Frame %d has invalid frame index in %s. Skipping.', frameIdx, csvFiles(i).name);
                continue;
            end
            
            % Initialize a structure for the current frame
            currentFrame = struct();
            
            for bp = 1:numBodyParts
                % Calculate the starting index for this body part's data
                startIdx = 1 + (bp-1)*3 + 1; % +1 to skip frame index
                
                % Extract x, y, and likelihood
                x = numericData(startIdx);
                y = numericData(startIdx + 1);
                likelihood = numericData(startIdx + 2);
                
                % Assign to the structure
                currentFrame.(bodyPartIdentifiers{bp}) = struct(...
                    'x', x, ...
                    'y', y, ...
                    'likelihood', likelihood ...
                );
            end
            
            % Add the current frame's data to the cell array
            bodyCoordinates{frameIdx} = currentFrame;
        end
    end
    
    % Close the CSV file
    fclose(fid);
    
    % Display total number of frames read
    fprintf('Total frames read from CSV: %d\n', frameIdx);
    
    %% Step 5: Define Body Parts List and Assign Weights
    
    % Define the list of body parts to plot
    bodyPartsList = {'Head_1', 'MainBody_2', 'Butt_3', 'TailTip_4', ...
                    'LeftForehand_5', 'RightForehand_6', 'LeftHindpaw_7', 'RightHindpaw_8'};
    
    % Define Weights for Each Body Part
    % Head: 30%, MainBody: 20%, Butt: 10%, TailTip: 0%
    % Each Forehand and Hindpaw: 10%
    weights = struct();
    weights.Head_1 = 30;
    weights.MainBody_2 = 20;
    weights.Butt_3 = 10;
    weights.TailTip_4 = 0;
    weights.LeftForehand_5 = 10;
    weights.RightForehand_6 = 10;
    weights.LeftHindpaw_7 = 10;
    weights.RightHindpaw_8 = 10;
    
    % Verify that total weights sum to 100%
    totalWeights = 0;
    bodyPartsNames = fieldnames(weights);
    for j = 1:length(bodyPartsNames)
        totalWeights = totalWeights + weights.(bodyPartsNames{j});
    end
    if totalWeights ~= 100
        warning('Total weights do not sum to 100%% for %s. Please check the weights.', csvFiles(i).name);
    else
        fprintf('Total weights assigned to body parts: %d%%\n', totalWeights);
    end
    
    %% Step 6: Calculate Total Speeds for Each Body Part in Centimeters
    
    % Create a time vector in seconds
    time = (0:frameIdx-1) / fps;
    
    % Preallocate structures to store speeds
    speeds = struct();
    
    % Define Speed Threshold and Continuous Duration
    speedThreshold = 10; % cm/s (updated from 20 cm/s)
    continuousDuration = 1; % seconds
    N = ceil(continuousDuration * fps); % Number of consecutive frames
    
    %% Calculate Total Speeds for Each Body Part in Centimeters
    
    for bp = 1:length(bodyPartsList)
        bodyPart = bodyPartsList{bp};
        
        % Initialize arrays to store x and y positions in cm
        x_positions_cm = zeros(1, frameIdx);
        y_positions_cm = zeros(1, frameIdx);
        
        % Extract x and y positions for the current body part and convert to cm
        for frame = 1:frameIdx
            currentFrame = bodyCoordinates{frame};
            
            if isfield(currentFrame, bodyPart)
                x = currentFrame.(bodyPart).x;
                y = currentFrame.(bodyPart).y;
                
                % Handle invalid or missing data
                if isnan(x) || isnan(y)
                    x_positions_cm(frame) = NaN;
                    y_positions_cm(frame) = NaN;
                else
                    x_positions_cm(frame) = x * pixels_to_cm;
                    y_positions_cm(frame) = y * pixels_to_cm;
                end
            else
                % If the body part is missing in this frame
                x_positions_cm(frame) = NaN;
                y_positions_cm(frame) = NaN;
            end
        end
        
        % Calculate total speed as the Euclidean distance between consecutive positions multiplied by FPS
        % Total speed will have one less element than the number of frames
        speed_total = sqrt(diff(x_positions_cm).^2 + diff(y_positions_cm).^2) * fps;
        
        % Store total speeds in the structure
        speeds.(bodyPart).speed_total = speed_total;
        
        % Also, create a time vector for speeds (since diff reduces length by 1)
        speeds.(bodyPart).time = time(2:end);
    end
    
    %% Step 7: Detect Continuous Immobilization
    
    for bp = 1:length(bodyPartsList)
        bodyPart = bodyPartsList{bp};
        
        % Initialize immobilization flags
        immobilized = false(1, length(speeds.(bodyPart).speed_total));
        
        % Detect immobilization for total speed
        if all(~isnan(speeds.(bodyPart).speed_total))
            binary = speeds.(bodyPart).speed_total < speedThreshold;
            conv_result = conv(binary, ones(1, N), 'same');
            immobilized = conv_result >= N;
        else
            warning('NaN values found in speed_total for %s in %s. Immobilization flags may be inaccurate.', bodyPart, csvFiles(i).name);
        end
        
        % Store immobilization flags
        speeds.(bodyPart).immobilized = immobilized;
    end
    
    %% Step 8: Plotting and Saving Results
    
    % Create a directory to save results if it doesn't exist
    resultsDir = fullfile(csvFiles(i).folder, 'Results');
    if ~exist(resultsDir, 'dir')
        mkdir(resultsDir);
    end
    
    %% 8.1: Plot and Save Total Speeds for All Body Parts as TIFF
    
    figure('Name', ['Mouse ' num2str(i) ' - Total Speeds for All Body Parts'], ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 1200, 600]); % 2:1 aspect ratio (width:height)
    
    for bp = 1:length(bodyPartsList)
        bodyPart = bodyPartsList{bp};
        
        % Select subplot position
        subplot(length(bodyPartsList),1,bp);
        
        % Plot Total Speed vs Time
        plot(speeds.(bodyPart).time, speeds.(bodyPart).speed_total, 'b', 'LineWidth', 1.5);
        hold on;
        yline(speedThreshold, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Threshold (%d cm/s)', speedThreshold));
        
        % Highlight immobilized periods
        plot(speeds.(bodyPart).time(speeds.(bodyPart).immobilized), ...
             speeds.(bodyPart).speed_total(speeds.(bodyPart).immobilized), ...
             'ro', 'MarkerSize', 4, 'DisplayName', 'Immobilized');
        
        xlabel('Time (s)');
        ylabel('Total Speed (cm/s)');
        title(sprintf('%s - Total Speed', bodyPart));
        legend('Total Speed', sprintf('Threshold (%d cm/s)', speedThreshold), 'Immobilized', 'Location', 'best');
        grid on;
        hold off;
    end
    
    % Adjust layout for better visibility
    tightfig();
    
    % Save the figure as TIFF
    [~, baseName, ~] = fileparts(csvFiles(i).name);
    totalSpeedsFileName = fullfile(resultsDir, [baseName '_Total_Speeds.tiff']);
    saveas(gcf, totalSpeedsFileName);
    fprintf('Total speeds plot saved to %s\n', totalSpeedsFileName);
    close;
    
    %% 8.2: Generate and Save Heatmap as TIFF
    
    % Prepare data for heatmap
    % Rows: Body Parts
    % Columns: Time Points
    speedMatrix = NaN(length(bodyPartsList), length(time)-1);
    
    for bp = 1:length(bodyPartsList)
        bodyPart = bodyPartsList{bp};
        speedMatrix(bp, :) = speeds.(bodyPart).speed_total;
    end
    
    % Create heatmap figure with 2:1 aspect ratio
    figure('Name', ['Mouse ' num2str(i) ' - Heatmap of Total Speeds'], ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 1200, 600]); % 2:1 aspect ratio
    
    % Use imagesc for heatmap
    imagesc(time(2:end), 1:length(bodyPartsList), speedMatrix);
    colorbar;
    colormap('hot');
    % Set color axis limits (adjust as needed)
    caxis([0 100]); 
    % Set Y-axis labels to body parts
    set(gca, 'YTick', 1:length(bodyPartsList), 'YTickLabel', bodyPartsList);
    
    xlabel('Time (s)');
    ylabel('Body Parts');
    title('Heatmap of Total Speeds (cm/s) for All Body Parts');
    grid on;
    
    % Overlay threshold lines on heatmap
    hold on;
    for bp = 1:length(bodyPartsList)
        plot([time(2) time(end)], [bp bp], 'w--', 'LineWidth', 0.5);
    end
    hold off;
    
    % Save the heatmap as TIFF
    heatmapFileName = fullfile(resultsDir, [baseName '_Speed_Heatmap.tiff']);
    saveas(gcf, heatmapFileName);
    fprintf('Heatmap saved to %s\n', heatmapFileName);
    close;
    
    %% 8.3: Calculate Weighted Immobility Score
    
    % Initialize an array to store total immobility scores per frame
    immobilityScore = zeros(1, frameIdx -1); % Since speed_total has one less frame due to diff
    
    % Iterate over each frame to calculate the weighted immobility score
    for frame = 1:(frameIdx -1)
        totalScore = 0;
        for bp = 1:length(bodyPartsList)
            bodyPart = bodyPartsList{bp};
            if speeds.(bodyPart).immobilized(frame)
                totalScore = totalScore + weights.(bodyPart);
            end
        end
        immobilityScore(frame) = totalScore;
    end
    
    % Add immobilityScore to the speeds structure
    speeds.immobilityScore = immobilityScore;
    
    % Define a threshold for immobility (e.g., 70%)
    immobilityThreshold = 70;
    immobileFrames = immobilityScore >= immobilityThreshold;
    
    %% 8.4: Plot and Save Weighted Immobility Score Over Time as TIFF
    
    figure('Name', ['Mouse ' num2str(i) ' - Weighted Immobility Score'], ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 1200, 600]); % 2:1 aspect ratio
    
    plot(time(2:end), immobilityScore, 'b', 'LineWidth', 1.5);
    hold on;
    yline(immobilityThreshold, 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf('Threshold (%d%%)', immobilityThreshold));
    
    % Highlight immobilized periods
    immobileTimes = time(2:end);
    plot(immobileTimes(immobileFrames), immobilityScore(immobileFrames), ...
        'ro', 'MarkerSize', 4, 'DisplayName', 'Immobilized');
    
    xlabel('Time (s)');
    ylabel('Immobility Score (%)');
    title('Total Weighted Immobility Score Over Time');
    legend('Immobility Score', sprintf('Threshold (%d%%)', immobilityThreshold), 'Immobilized', 'Location', 'best');
    grid on;
    
    % Set x-axis upper limit to 360 seconds
    xlim([0, 360]);
    
    hold off;
    
    % Calculate summary statistics for immobility
    totalImmobilizedTime = sum(immobileFrames) / fps;
    totalTime = frameIdx / fps;
    percentageImmobilized = (sum(immobileFrames) / (frameIdx -1)) * 100;
    
    fprintf('Total Immobilized Time: %.2f seconds out of %.2f seconds (%.2f%%)\n', ...
        totalImmobilizedTime, totalTime, percentageImmobilized);
    
    % Optionally, add a text box to the immobility score plot
    dim = [.15 .6 .3 .3];
    str = sprintf('Total Immobilized Time: %.2f s (%.2f%%)', totalImmobilizedTime, percentageImmobilized);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'white');
    
    % Save the immobility score plot as TIFF
    immobilityScoreFileName = fullfile(resultsDir, [baseName '_Immobility_Score.tiff']);
    saveas(gcf, immobilityScoreFileName);
    fprintf('Immobility score plot saved to %s\n', immobilityScoreFileName);
    close;
    
    %% 8.5: Generate and Save Annotated Output Video
    
    fprintf('Generating annotated video for Mouse %d...\n', i);
    
    % Create a VideoReader object for the scaling MP4
    vidReader = VideoReader(videoFileFullPath);
    
    % Define the output video file name
    annotatedVideoName = fullfile(resultsDir, [baseName '_Annotated.mp4']);
    
    % Create a VideoWriter object
    vidWriter = VideoWriter(annotatedVideoName, 'MPEG-4');
    vidWriter.FrameRate = fps;
    open(vidWriter);
    
    % Reset the VideoReader to the beginning
    vidReader.CurrentTime = 0;
    
    % Determine the number of frames to process (up to 360 seconds)
    maxFrames = min(ceil(360 * fps), frameIdx);
    
    % Loop through each frame and add annotations
    for frameNumber = 1:maxFrames
        if hasFrame(vidReader)
            frame = readFrame(vidReader);
        else
            % If the video has fewer frames than expected, pad with the last frame
            frame = lastFrame;
        end
        
        lastFrame = frame; % Store the last frame in case video ends early
        
        % Determine the immobility status for the current frame
        % Note: immobilityScore has (frameIdx -1) frames, so map accordingly
        if frameNumber <= length(immobileFrames)
            isImmobile = immobileFrames(frameNumber);
        else
            isImmobile = false; % Default to moving if data is insufficient
        end
        
        % Define annotation properties
        if isImmobile
            % Red dot and "immobile" text
            circleColor = 'red';
            statusText = 'immobile';
        else
            % Green dot and "moving" text
            circleColor = 'green';
            statusText = 'moving';
        end
        
        % Define the position for the symbol and text (top right corner)
        [frameHeight, frameWidth, ~] = size(frame);
        margin = 20; % pixels from the edge
        circleRadius = 10; % pixels
        circlePosition = [frameWidth - margin - circleRadius, margin + circleRadius, circleRadius];
        textPosition = [frameWidth - margin - 150, margin - 10]; % Adjust as needed
        
        % Insert the colored circle
        frame = insertShape(frame, 'FilledCircle', circlePosition, 'Color', circleColor, 'Opacity', 1);
        
        % Insert the status text
        frame = insertText(frame, textPosition, statusText, 'FontSize', 18, 'BoxColor', 'white', ...
                          'BoxOpacity', 0.6, 'TextColor', 'black', 'AnchorPoint', 'LeftTop');
        
        % Write the annotated frame to the output video
        writeVideo(vidWriter, frame);
    end
    
    % Close the VideoWriter
    close(vidWriter);
    fprintf('Annotated video saved to %s\n', annotatedVideoName);
    
end

%% Step 10: Save Summary Statistics to Excel

if ~isempty(summaryData)
    % Convert summaryData to table
    summaryTable = cell2table(summaryData, 'VariableNames', {'Mouse', 'Total_Immobilized_Time_s', 'Latency_to_First_Immobility_s'});
    
    % Define the Excel file name
    excelFileName = 'Summary_TST.xlsx'; % Changed from FST to TST
    
    % Write the table to Excel
    try
        writetable(summaryTable, excelFileName, 'Sheet', 1);
        fprintf('\nSummary statistics saved to %s\n', excelFileName);
    catch ME
        warning('Failed to write summary statistics to Excel: %s', ME.message);
    end
else
    warning('No summary data to save.');
end

%% Helper Function: Adjusts figure layout to minimize white space

function tightfig()
    % Adjust the figure to remove excess white space
    set(gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
    ax = findall(gcf, 'type', 'axes');
    for i = 1:length(ax)
        outerpos = ax(i).OuterPosition;
        ti = ax(i).TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax(i).Position = [left bottom ax_width ax_height];
    end
end