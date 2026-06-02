% MATLAB Script to Overlay Body Part Coordinates on Video Frames and Save Coordinates 

% Clear workspace, close all figures, and clear command window
close all;
clear; clc;

%% Parameters
% Define file names (adjust these as needed)
csvFileName = 'coordinates.csv';                % Your CSV file with coordinates
videoFileName = 'input_video.mp4';             % Your original video file
outputVideoName = 'output_with_overlays.mp4'; % Output video with overlays
coordinatesSaveFile = 'bodyCoordinates.mat';    % File to save coordinates for future use

% Define visualization parameters
markerSize = 10;                               % Size of the markers
likelihoodThreshold = 0.05;                    % Minimum likelihood to display a marker

% Define whether to enable interactive frame-by-frame inspection
enableInteractive = true;

%% Step 1: Read and Parse the CSV Data

% Open the CSV file
fid = fopen(csvFileName, 'r');
if fid == -1
    error('Cannot open the CSV file: %s', csvFileName);
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
    error('Mismatch in the number of columns among the header lines.');
end

% Extract body part names and coordinate types (excluding the first column which is 'scorer')
bodyParts = bodyPartsHeaders(2:end); % Cell array (e.g., {'Head', 'Head', 'Head', 'MainBody', ...})
coordTypes = coordsHeaders(2:end);    % Should be {'x', 'y', 'likelihood', 'x', 'y', 'likelihood', ...}

% Determine the number of body parts
% Each body part has three columns: x, y, likelihood
if mod((numColumns -1), 3) ~= 0
    error('Each body part should have exactly three columns: x, y, likelihood.');
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
            warning('Frame %d has invalid frame index. Skipping.', frameIdx);
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

%% Step 2: Load the Video and Extract FPS

vidReader = VideoReader(videoFileName);

% Extract and round up the FPS
fps_original = vidReader.FrameRate;
fps = ceil(fps_original);
fprintf('Original FPS: %.2f\n', fps_original);
fprintf('Rounded FPS: %d\n', fps);

% Calculate the total number of video frames based on rounded FPS
totalVideoFrames = floor(vidReader.Duration * fps);
fprintf('Total frames in video: %d\n', totalVideoFrames);

% Check if CSV has more frames than the video
if frameIdx > totalVideoFrames
    warning('CSV has more frames (%d) than the video (%d). Truncating to video length.', frameIdx, totalVideoFrames);
    bodyCoordinates = bodyCoordinates(1:totalVideoFrames);
    numFrames = totalVideoFrames;
else
    numFrames = frameIdx;
end

%% Step 3: User Input for Start and End Times

% Prompt user for start and end times in seconds
fprintf('Video Duration: %.2f seconds\n', vidReader.Duration);

validInput = false;
while ~validInput
    startTime = input('Enter the start time in seconds (>=0): ');
    endTime = input('Enter the end time in seconds (<= video duration): ');
    
    % Validate inputs
    if isnumeric(startTime) && isnumeric(endTime) && ...
       startTime >= 0 && endTime <= vidReader.Duration && ...
       startTime < endTime
        validInput = true;
    else
        disp('Invalid input. Please ensure that 0 <= start time < end time <= video duration.');
    end
end

% Convert start and end times to frame numbers based on rounded FPS
startFrame = ceil(startTime * fps);
endFrame = ceil(endTime * fps);

% Ensure frame numbers are within bounds
startFrame = max(startFrame, 1);
endFrame = min(endFrame, numFrames);

fprintf('Processing frames from %d to %d.\n', startFrame, endFrame);

% Update the number of frames to process
numFramesToProcess = endFrame - startFrame + 1;

%% Step 4: Save Coordinates for Future Use

% Save the bodyCoordinates cell array to a .mat file
save(coordinatesSaveFile, 'bodyCoordinates');
fprintf('Body coordinates saved to %s\n', coordinatesSaveFile);

%% Step 5: Initialize Video Writer

vidWriter = VideoWriter(outputVideoName, 'MPEG-4');
vidWriter.FrameRate = fps;
open(vidWriter);

%% Step 6: Skip Frames to Reach Start Frame

if startFrame > 1
    fprintf('Skipping the first %d frames to reach the start frame.\n', startFrame-1);
    for skipIdx = 1:startFrame-1
        if hasFrame(vidReader)
            readFrame(vidReader); % Read and discard frames
        else
            error('Cannot skip to startFrame: video has fewer frames.');
        end
    end
end

%% Step 7: Process and Overlay Coordinates

% Create a figure for interactive inspection if enabled
if enableInteractive
    hFig = figure('Name', 'Frame Inspection', 'NumberTitle', 'off');
end

% Define a single color for all body parts (e.g., green)
defaultColor = [0, 1, 0]; % RGB for green

for idx = 1:numFramesToProcess
    currentFrameIdx = startFrame + idx - 1;
    
    % Read the current frame from the video
    if hasFrame(vidReader)
        frame = readFrame(vidReader);
    else
        warning('Reached end of video before processing all frames.');
        break;
    end
    
    % Get coordinates for the current frame
    currentFrame = bodyCoordinates{currentFrameIdx};
    
    % Overlay coordinates on the frame
    overlaidFrame = frame;
    
    for bp = 1:numBodyParts
        bodyPartID = bodyPartIdentifiers{bp};
        x = currentFrame.(bodyPartID).x;
        y = currentFrame.(bodyPartID).y;
        likelihood = currentFrame.(bodyPartID).likelihood;
        
        % Only plot if likelihood is above the threshold
        if likelihood >= likelihoodThreshold
            % Extract base body part name by removing the suffix (e.g., 'Head_1' -> 'Head')
            underscoreIdx = strfind(bodyPartID, '_');
            if ~isempty(underscoreIdx)
                baseBodyPart = bodyPartID(1:underscoreIdx(1)-1);
            else
                baseBodyPart = bodyPartID;
            end
            
            % Draw a filled circle at (x, y) with the default color
            overlaidFrame = insertShape(overlaidFrame, 'FilledCircle', [x, y, markerSize], ...
                'Color', defaultColor, 'Opacity', 1);
            
            % Add text label near the marker
            overlaidFrame = insertText(overlaidFrame, [x + 5, y + 5], baseBodyPart, ...
                'TextColor', defaultColor, 'FontSize', 8, 'BoxOpacity', 0);
        end
    end
    
    % Write the overlaid frame to the output video
    writeVideo(vidWriter, overlaidFrame);
    
    % If interactive mode is enabled, display the frame for inspection
    if enableInteractive
        imshow(overlaidFrame);
        title(sprintf('Frame %d / %d', currentFrameIdx, endFrame));
        drawnow;
        % Wait for user to press a key to continue
        pause(0.01); % Brief pause to allow figure to update
    end
    
    % Optional: Display progress every 100 frames
    if mod(idx, 100) == 0
        fprintf('Processed %d/%d frames.\n', idx, numFramesToProcess);
    end
end

% Close the video writer
close(vidWriter);
disp(['Overlaid video saved as ', outputVideoName]);

%% Optional: Close the interactive figure
if enableInteractive
    close(hFig);
end
