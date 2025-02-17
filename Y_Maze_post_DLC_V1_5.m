clear all; 
close all;
clc;

%% Set up parameters
% Likelihood threshold (if filtering is still needed)
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

% Analyze first 8 minutes of data
analysis_duration_sec = 8 * 60;  % 480 seconds
max_frames_to_analyze = frames_per_second * analysis_duration_sec;

%% --- Define Y-Maze Arms ROIs ---
% Use the first frame of the first video for ROI selection
firstFrame = readFrame(firstVideoObj);
figure('Name','Draw Y-maze Arms ROIs');
imshow(firstFrame);

% Draw polygon for A-arm
title('Draw polygon for A-arm, then double-click to finish');
roi_A = drawpolygon;
wait(roi_A);
mask_A = roi_A.Position;

% Draw polygon for B-arm
title('Draw polygon for B-arm, then double-click to finish');
roi_B = drawpolygon;
wait(roi_B);
mask_B = roi_B.Position;

% Draw polygon for C-arm
title('Draw polygon for C-arm, then double-click to finish');
roi_C = drawpolygon;
wait(roi_C);
mask_C = roi_C.Position;
close(gcf);

%% --- Process Each File ---
% Initialize summary data: CSV file, arm sequence, and alternation percentage
summaryData = {};

for i = 1:length(csv_files)
    csv_fileName = csv_files(i).name;
    video_fileName = video_files(i).name;
    fprintf('Processing file: %s\n', csv_fileName);

    % Read CSV data
    csvData = readtable(csv_fileName);
    % Assuming the CSV contains:
    %   Column 1: Frame numbers
    %   Column 8: MainBody X coordinates
    %   Column 9: MainBody Y coordinates
    %   Column 10: Likelihood values
    frames = csvData{:, 1};
    MainBodyX = csvData{:, 8};
    MainBodyY = csvData{:, 9};
    likelihood = csvData{:, 10};

    % Limit analysis to first 8 minutes or available frames
    num_frames = min(max_frames_to_analyze, length(frames));
    frames = frames(1:num_frames);
    MainBodyX = MainBodyX(1:num_frames);
    MainBodyY = MainBodyY(1:num_frames);
    likelihood = likelihood(1:num_frames);

    % Initialize sequence and state tracking
    sequence = '';   % To record the sequence (e.g., 'ABCCBA')
    current_state = 'O';  % 'O' stands for "outside" any ROI

    for j = 1:num_frames
        % Use detection only if likelihood is high; otherwise treat as 'O'
        if likelihood(j) < likelihood_threshold
            new_state = 'O';
        else
            % Get current coordinate
            x = MainBodyX(j);
            y = MainBodyY(j);

            % Determine which arm the point falls into
            inA = inpolygon(x, y, mask_A(:,1), mask_A(:,2));
            inB = inpolygon(x, y, mask_B(:,1), mask_B(:,2));
            inC = inpolygon(x, y, mask_C(:,1), mask_C(:,2));

            if inA
                new_state = 'A';
            elseif inB
                new_state = 'B';
            elseif inC
                new_state = 'C';
            else
                new_state = 'O';
            end
        end

        % If the state changes, and the new state is an arm (not 'O'),
        % record it as a new entry.
        if new_state ~= current_state
            if new_state ~= 'O'
                sequence = [sequence new_state];  %#ok<AGROW>
            end
            current_state = new_state;
        end
    end

    %% --- Calculate Spontaneous Alternation Percentage ---
    % An alternation is defined as any set of three consecutive entries that
    % contain all three arms (A, B, C).
    alternation_count = 0;
    total_triplets = length(sequence) - 2;

    if total_triplets > 0
        for j = 1:total_triplets
            triplet = sequence(j:j+2);
            if numel(unique(triplet)) == 3
                alternation_count = alternation_count + 1;
            end
        end
        alternation_percentage = (alternation_count / total_triplets) * 100;
    else
        alternation_percentage = NaN; % Not enough entries to compute
    end

    %% --- Plot Trajectory on First Frame with ROIs ---
    % Reload the first frame from the current video
    videoObj = VideoReader(video_fileName);
    firstFrameCurrent = readFrame(videoObj);
    figure;
    imshow(firstFrameCurrent);
    hold on;

    % Plot the ROI outlines for A, B, and C arms
    % Close the polygon by repeating the first point at the end.
    plot([mask_A(:,1); mask_A(1,1)], [mask_A(:,2); mask_A(1,2)], 'r-', 'LineWidth', 2, 'DisplayName', 'A-arm');
    plot([mask_B(:,1); mask_B(1,1)], [mask_B(:,2); mask_B(1,2)], 'b-', 'LineWidth', 2, 'DisplayName', 'B-arm');
    plot([mask_C(:,1); mask_C(1,1)], [mask_C(:,2); mask_C(1,2)], 'y-', 'LineWidth', 2, 'DisplayName', 'C-arm');

    % Extract valid trajectory points
    validPoints = likelihood >= likelihood_threshold;
    x_traj = MainBodyX(validPoints);
    y_traj = MainBodyY(validPoints);

    % Plot continuous trajectory line
    plot(x_traj, y_traj, 'g-', 'LineWidth', 2, 'DisplayName', 'Trajectory');

    % Add arrowheads to indicate direction
    numArrows = 20;  % Adjust the number of arrows as needed
    totalPoints = length(x_traj);
    if totalPoints > 1
        indices = round(linspace(1, totalPoints-1, numArrows));
        dx = diff(x_traj);
        dy = diff(y_traj);
        quiver(x_traj(indices), y_traj(indices), dx(indices), dy(indices), 0, ...
            'MaxHeadSize', 2, 'Color', 'g', 'LineWidth', 1, 'AutoScale', 'off');
    end

    title(sprintf('Video: %s - Trajectory with Y-maze Arms', video_fileName), 'Interpreter', 'none');
    legend('show','Location','best');
    hold off;

    % Save the trajectory image as a TIFF file
    [~, videoFileNameNoExt, ~] = fileparts(video_fileName);
    saveas(gcf, sprintf('%s_Trajectory_FirstFrame.tif', videoFileNameNoExt));
    close(gcf);

    %% --- Append Results to Summary ---
    summaryData(end+1,:) = {csv_fileName, sequence, alternation_percentage}; %#ok<SAGROW>

    fprintf('Processed %s: Sequence: %s, Alternation%%: %.2f\n', csv_fileName, sequence, alternation_percentage);
end

%% --- Save Summary Data ---
summaryTable = cell2table(summaryData, 'VariableNames', {'CSV_File','Sequence','Alternation_Percentage'});
writetable(summaryTable, 'Summary_Y_maze.xlsx');
fprintf('Summary data saved to Summary_Y_maze.xlsx\n');
