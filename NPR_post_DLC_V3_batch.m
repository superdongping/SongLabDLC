clear all;  
close all;
clc;

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

% Let the user draw two circles for the objects inside the open field
title('Now, draw the two circular ROIs around the objects');
roi1 = drawcircle('Color','r'); % User draws the first circle
roi2 = drawcircle('Color','r'); % User draws the second circle, same radius

% Get the center and radius of the drawn circles
center1 = roi1.Center;
radius1 = roi1.Radius;
center2 = roi2.Center;
radius2 = radius1; % Ensure both circles have the same radius

% Parameters for exploration detection
window_duration = 0.2; % Seconds (changed from 0.5)
min_overlap = 0.8; % Minimum fraction of frames within ROI to consider a valid window
movement_threshold = 1.5; % Pixels (lowered from 2)
stable_frames_threshold = 2 * 30; % 2 seconds, assuming default fps of 30

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
    
    % Calculate distances from nose to ROI centers
    dist1 = sqrt((noseX - center1(1)).^2 + (noseY - center1(2)).^2);
    dist2 = sqrt((noseX - center2(1)).^2 + (noseY - center2(2)).^2);
    
    % Determine if nose is within ROI
    inROI1 = dist1 <= radius1;
    inROI2 = dist2 <= radius2;
    
    % Detect exploration for both ROIs
    [exploreROI1, exploreEventsROI1, totalTimeROI1, timestampsROI1] = detectExploration(inROI1, center1, radius1, bodyX, bodyY, noseX, noseY, fps, window_size, min_overlap, movement_threshold);
    [exploreROI2, exploreEventsROI2, totalTimeROI2, timestampsROI2] = detectExploration(inROI2, center2, radius2, bodyX, bodyY, noseX, noseY, fps, window_size, min_overlap, movement_threshold);
    
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

% Save the summary data to an Excel file
summaryTable = cell2table(summaryData, 'VariableNames', {'CSV_File', 'Time_Object_1_sec', 'Events_Object_1', 'Time_Object_2_sec', 'Events_Object_2'});
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

function [exploreFlags, numEvents, totalTime, timestamps] = detectExploration(inROI, center, radius, bodyX, bodyY, noseX, noseY, fps, window_size, min_overlap, movement_threshold)
    % Detect exploration events within an ROI using a sliding window approach
    %
    % Parameters:
    %   inROI            - Logical array indicating if the nose is within the ROI
    %   center           - [x, y] coordinates of the ROI center
    %   radius           - Radius of the ROI
    %   bodyX, bodyY     - Coordinates of the mouse's body
    %   noseX, noseY     - Coordinates of the mouse's nose
    %   fps              - Frames per second of the video
    %   window_size      - Number of frames in the exploration window
    %   min_overlap      - Minimum fraction of frames within ROI to consider a valid window
    %   movement_threshold - Maximum allowed standard deviation in position to consider as stable
    
    % Initialize exploration flags
    exploreFlags = false(length(inROI),1);
    numEvents = 0;
    totalTime = 0;
    timestamps = []; % To store exploration event start times
    inEvent = false;
    
    % Calculate the minimum number of frames within ROI required
    min_frames_in_ROI = ceil(min_overlap * window_size);
    
    % Loop through frames with a sliding window
    j = 1;
    while j <= (length(inROI) - window_size +1)
        windowEnd = j + window_size -1;
        current_window = inROI(j:windowEnd);
        num_in_ROI = sum(current_window);
        
        if num_in_ROI >= min_frames_in_ROI
            % Calculate distance at start and end of window
            distance_start = sqrt((noseX(j) - center(1)).^2 + (noseY(j) - center(2)).^2);
            distance_end = sqrt((noseX(windowEnd) - center(1)).^2 + (noseY(windowEnd) - center(2)).^2);
            
            % Check if distance has decreased
            if distance_end < distance_start
                % Potential exploration event detected at frame j
                if ~inEvent
                    % Check if not sitting on the object
                    body_dist = sqrt((bodyX(j) - center(1)).^2 + (bodyY(j) - center(2)).^2);
                    overlap = body_dist <= (radius / 2);
                    
                    if overlap
                        % Check if the body is stable for a certain duration
                        stable = checkStability(bodyX, bodyY, j, stable_frames_threshold, movement_threshold);
                    else
                        stable = false;
                    end
                    
                    if ~stable
                        % Mark exploration event
                        exploreFlags(j:windowEnd) = true;
                        numEvents = numEvents +1;
                        totalTime = totalTime + window_size / fps;
                        timestamps = [timestamps; (j / fps)]; % Convert frame to time in seconds
                        inEvent = true;
                        % Skip the window to avoid overlapping events
                        j = windowEnd;
                        continue; % Continue to next iteration
                    end
                end
            end
        end
        
        % Reset inEvent flag if not in an event
        if exploreFlags(j) == false
            inEvent = false;
        end
        
        % Move to the next frame
        j = j +1;
    end
end

function isStable = checkStability(bodyX, bodyY, startIdx, window, movement_threshold)
    % Check if the body position is stable over a window of frames
    %
    % Parameters:
    %   bodyX              - X coordinates of the body
    %   bodyY              - Y coordinates of the body
    %   startIdx           - Starting index of the stability check
    %   window             - Number of frames to check for stability
    %   movement_threshold - Maximum allowed standard deviation in position
    
    if startIdx + window -1 > length(bodyX)
        isStable = false;
        return;
    end
    
    windowX = bodyX(startIdx:startIdx+window-1);
    windowY = bodyY(startIdx:startIdx+window-1);
    
    if std(windowX) < movement_threshold && std(windowY) < movement_threshold
        isStable = true;
    else
        isStable = false;
    end
end
