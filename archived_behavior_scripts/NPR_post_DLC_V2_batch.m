clear all; 
close all;
clc;

% Set the likelihood threshold
likelihood_threshold = 0.05;

% Search for CSV and MP4 files in the current folder
csv_files = dir('*.csv');
video_files = dir('*.mp4');

if isempty(csv_files) || isempty(video_files)
    error('No CSV or MP4 files found in the current directory.');
end

% Load the first video to define the open field and circles
video_fileName = video_files(1).name; 
fprintf('Loading first video: %s\n', video_fileName);

% Load the first video and get the first frame
videoObj = VideoReader(video_fileName);
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

% Loop through all videos and CSV files
summaryData = [];
for i = 1:length(csv_files)
    % Load the CSV file
    csv_fileName = csv_files(i).name;
    fprintf('Reading CSV file: %s\n', csv_fileName);
    csvData = readtable(csv_fileName);
    
    % Load the video file
    video_fileName = video_files(i).name;
    fprintf('Loading video file: %s\n', video_fileName);
    videoObj = VideoReader(video_fileName);
    
    % Get the frame information from the CSV
    frames = csvData{:, 1};      % 1st column: Frame numbers
    noseX = csvData{:, 2};       % 2nd column: X coordinates of the head
    noseY = csvData{:, 3};       % 3rd column: Y coordinates of the head
    likelihood = csvData{:, 4};  % 4th column: Likelihood/confidence of the head detection
    
    % Exclude points with low likelihood
    validPoints = likelihood >= likelihood_threshold;
    noseX = noseX(validPoints);
    noseY = noseY(validPoints);
    
    % Exclude points outside the open field rectangle
    inOpenField = noseX >= openFieldPosition(1) & noseX <= (openFieldPosition(1) + openFieldPosition(3)) & ...
                  noseY >= openFieldPosition(2) & noseY <= (openFieldPosition(2) + openFieldPosition(4));
    noseX = noseX(inOpenField);
    noseY = noseY(inOpenField);
    
    % Plot the trajectory on the first frame of this video
    firstFrame = readFrame(videoObj);
    figure;
    imshow(firstFrame);
    hold on;
    plot(noseX, noseY, 'm', 'LineWidth', 1); % Plot the nose trajectory
    viscircles(center1, radius1, 'EdgeColor', 'r'); % Draw the first circle
    viscircles(center2, radius2, 'EdgeColor', 'r'); % Draw the second circle
    rectangle('Position', openFieldPosition, 'EdgeColor', 'g', 'LineWidth', 2); % Draw the open field rectangle
    title(sprintf('Video: %s', video_fileName));
    hold off;
    
    % Save the frame with circles, trajectory, and open field as a TIFF file
    [~, videoFileNameNoExt, ~] = fileparts(video_fileName); % Remove extension from video file name
    saveas(gcf, sprintf('%s_FirstFrame_with_ROIs.tif', videoFileNameNoExt));

    % Calculate the distance from the nose to each object
    dist1 = sqrt((noseX - center1(1)).^2 + (noseY - center1(2)).^2);
    dist2 = sqrt((noseX - center2(1)).^2 + (noseY - center2(2)).^2);
    
    % Count the number of frames where the nose is within each ROI
    framesInROI1 = sum(dist1 <= radius1);
    framesInROI2 = sum(dist2 <= radius2);
    
    % Calculate the total time spent investigating each object
    fps = 30; % Assuming 30 frames per second
    timeInROI1 = framesInROI1 / fps;
    timeInROI2 = framesInROI2 / fps;
    
    % Append the results to the summary data
    summaryData = [summaryData; {csv_fileName, timeInROI1, timeInROI2}];
    
    % Display the results for this video
    fprintf('Video %d - Time spent on object 1: %.2f seconds\n', i, timeInROI1);
    fprintf('Video %d - Time spent on object 2: %.2f seconds\n', i, timeInROI2);
end

% Save the summary data to an Excel file
summaryTable = cell2table(summaryData, 'VariableNames', {'CSV_File', 'Time_Object_1', 'Time_Object_2'});
writetable(summaryTable, 'summary_NPR.xlsx');
fprintf('Summary data saved to summary_NPR.xlsx\n');
