% fiber_photometry_behavior_analysis_combined_sorted_modified.m

% This script analyzes and visualizes multiple behavior-tagged fiber photometry datasets.
% For each CSV and corresponding Excel file pair, it performs the following:
%   1. Imports raw CSV data and behavior tags from an Excel file.
%   2. Applies a low-pass Butterworth filter to smooth the Real_Signal.
%   3. Aligns the smoothed signal to behavior onset times.
%   4. Computes ΔF/F normalization.
%   5. Collects normalized data for consistent plotting.
%   6. After processing all file pairs, it generates:
%       - Individual overlay plots with mean and SEM using consistent y-axis ranges.
%       - Individual heatmaps with unified color scales.
%       - Aggregate overlay plots and sorted heatmaps based on AUC from 0 to 100s post-tag.

function fiber_photometry_behavior_analysis_combined_sorted_modified()
    % Clear workspace and close all figures
    close all;
    clearvars;
    warning off;

    %% 1. Import Multiple CSV and Excel Files

    % Prompt user to select multiple CSV raw data files
    [csvFileNames, csvPath] = uigetfile('*.csv', 'Select the Raw CSV Data Files', 'MultiSelect', 'on');
    if isequal(csvFileNames,0)
        disp('User canceled CSV file selection.');
        return;
    end

    % If only one file is selected, uigetfile returns a char, convert to cell array
    if ischar(csvFileNames)
        csvFileNames = {csvFileNames};
    end

    numCsvFiles = length(csvFileNames);
    fprintf('Number of CSV files selected: %d\n', numCsvFiles);

    % Prompt user to select multiple Excel Behavior Tags files
    [xlsxFileNames, xlsxPath] = uigetfile('*.xlsx', 'Select the Behavior Tags Excel Files', 'MultiSelect', 'on');
    if isequal(xlsxFileNames,0)
        disp('User canceled Excel file selection.');
        return;
    end

    % If only one file is selected, uigetfile returns a char, convert to cell array
    if ischar(xlsxFileNames)
        xlsxFileNames = {xlsxFileNames};
    end

    numXlsxFiles = length(xlsxFileNames);
    fprintf('Number of Excel files selected: %d\n', numXlsxFiles);

    % Check if number of CSV and Excel files match
    if numCsvFiles ~= numXlsxFiles
        error('Number of CSV files (%d) and Excel files (%d) must match.', numCsvFiles, numXlsxFiles);
    end

    %% 2. Initialize Aggregated Data Containers

    all_normalizedData1 = []; % To store normalized data for Behavior1 across all files
    all_normalizedData2 = []; % To store normalized data for Behavior2 across all files

    all_timeWindow1 = []; % To store time windows for Behavior1 (assuming same across files)
    all_timeWindow2 = []; % To store time windows for Behavior2 (assuming same across files)

    % Initialize global min and max variables for y-axis and heatmap color limits
    global_min1 = Inf;
    global_max1 = -Inf;
    global_min2 = Inf;
    global_max2 = -Inf;

    global_heatmap_min1 = Inf;
    global_heatmap_max1 = -Inf;
    global_heatmap_min2 = Inf;
    global_heatmap_max2 = -Inf;

    %% 3. Initialize Data Structure to Store Normalized Data Per File

    % Preallocate structure array for storing per-file normalized data
    fileData(numCsvFiles) = struct(...
        'csvFileName', '', ...
        'xlsxFileName', '', ...
        'normalizedData1', [], ...
        'timeWindow1', [], ...
        'normalizedData2', [], ...
        'timeWindow2', []);

    %% 4. First Pass: Process Each CSV and Excel File Pair Individually (Data Collection Only)

    for i = 1:numCsvFiles
        fprintf('\nProcessing File Pair %d of %d\n', i, numCsvFiles);

        % Get CSV and Excel file names
        csvFileName = csvFileNames{i};
        xlsxFileName = xlsxFileNames{i};

        % Display processing information
        fprintf('CSV File: %s\n', csvFileName);
        fprintf('Excel File: %s\n', xlsxFileName);

        % Full paths
        csvFullPath = fullfile(csvPath, csvFileName);
        xlsxFullPath = fullfile(xlsxPath, xlsxFileName);

        %% 5. Load and Process CSV Data

        try
            rawData = readtable(csvFullPath);
        catch ME
            warning('Failed to read CSV file %s: %s. Skipping this pair.', csvFileName, ME.message);
            continue;
        end

        % Check if required columns exist
        requiredColumns = [1, 2, 3];
        if width(rawData) < max(requiredColumns)
            warning('CSV file %s does not contain the required columns. Skipping this pair.', csvFileName);
            continue;
        end

        % Extract timestamps and LED signals
        try
            TimeStamp = rawData{:,1};  % Column 1
            LED_410 = rawData{:,2};    % Column 2
            LED_470 = rawData{:,3};    % Column 3
        catch ME
            warning('Error extracting data from CSV file %s: %s. Skipping this pair.', csvFileName, ME.message);
            continue;
        end

        % Validate and convert TimeStamp to seconds
        if isduration(TimeStamp)
            timeVector = seconds(TimeStamp) - seconds(TimeStamp(1));
        elseif isdatetime(TimeStamp)
            timeVector = seconds(TimeStamp - timeStamp(1));
        else
            warning('Timestamp column in CSV file %s must be of type datetime or duration. Skipping this pair.', csvFileName);
            continue;
        end
        timeVector = timeVector(:);  % Ensure column vector

        % Calculate Real_Signal (470/410)
        Real_Signal = LED_470 ./ LED_410;
        Real_Signal = Real_Signal(:);  % Ensure column vector

        % Handle non-finite values
        finiteIdx = isfinite(timeVector) & isfinite(Real_Signal);
        if any(~finiteIdx)
            warning('Non-finite values found in CSV file %s. These will be removed.', csvFileName);
            timeVector = timeVector(finiteIdx);
            Real_Signal = Real_Signal(finiteIdx);
        end

        % Check if data is sufficient after cleaning
        if isempty(timeVector) || isempty(Real_Signal)
            warning('No valid data available in CSV file %s after cleaning. Skipping this pair.', csvFileName);
            continue;
        end

        %% 6. Apply Low-Pass Butterworth Filter to Smooth the Signal

        % Define filter parameters
        Fc = 5;             % Cutoff frequency in Hz
        filterOrder = 4;    % Order of the Butterworth filter

        % Estimate the sampling frequency (Fs)
        samplingInterval = median(diff(timeVector));  % Median sampling interval in seconds
        Fs = 1 / samplingInterval;                     % Sampling frequency in Hz
        fprintf('Estimated Sampling Frequency: %.2f Hz\n', Fs);

        % Design Butterworth low-pass filter
        [b, a] = butter(filterOrder, Fc/(Fs/2), 'low');

        % Apply zero-phase filtering using filtfilt to prevent phase distortion
        Real_Signal_Smoothed = filtfilt(b, a, Real_Signal);

        %% 7. Load and Extract Behavior Tags from Excel File

        try
            [~, sheets] = xlsfinfo(xlsxFullPath);
            if isempty(sheets)
                warning('No sheets found in Excel file %s. Skipping behavior tag extraction.', xlsxFileName);
                behavior1Times = [];
                behavior2Times = [];
            else
                % Initialize behavior times
                behavior1Times = [];
                behavior2Times = [];

                % Read Behavior1 sheet if it exists
                if ismember('Behavior1', sheets)
                    try
                        behavior1Data = readtable(xlsxFullPath, 'Sheet', 'Behavior1');
                        if ismember('Behavior1_Onset_Time_s', behavior1Data.Properties.VariableNames)
                            behavior1Times = behavior1Data.Behavior1_Onset_Time_s;
                            fprintf('Loaded %d Behavior1 tags from %s.\n', length(behavior1Times), xlsxFileName);
                        else
                            warning('Sheet "Behavior1" in %s does not contain the column "Behavior1_Onset_Time_s".', xlsxFileName);
                        end
                    catch ME
                        warning('Failed to read "Behavior1" sheet in %s: %s', xlsxFileName, ME.message);
                    end
                else
                    warning('Sheet "Behavior1" not found in Excel file %s.', xlsxFileName);
                end

                % Read Behavior2 sheet if it exists
                if ismember('Behavior2', sheets)
                    try
                        behavior2Data = readtable(xlsxFullPath, 'Sheet', 'Behavior2');
                        if ismember('Behavior2_Onset_Time_s', behavior2Data.Properties.VariableNames)
                            behavior2Times = behavior2Data.Behavior2_Onset_Time_s;
                            fprintf('Loaded %d Behavior2 tags from %s.\n', length(behavior2Times), xlsxFileName);
                        else
                            warning('Sheet "Behavior2" in %s does not contain the column "Behavior2_Onset_Time_s".', xlsxFileName);
                        end
                    catch ME
                        warning('Failed to read "Behavior2" sheet in %s: %s', xlsxFileName, ME.message);
                    end
                else
                    warning('Sheet "Behavior2" not found in Excel file %s.', xlsxFileName);
                end
            end
        catch ME
            warning('Failed to read Excel file %s: %s. Skipping behavior tag extraction.', xlsxFileName, ME.message);
            behavior1Times = [];
            behavior2Times = [];
        end

        %% 8. Define Analysis Parameters

        % Define time window around behavior onset
        preTime = -30;  % seconds before tag
        postTime = 30; % seconds after tag
        windowDuration = postTime - preTime; % total window duration

        %% 9. Analyze Behavior1 Tags (First Pass)

        if ~isempty(behavior1Times)
            fprintf('Analyzing Behavior 1 tags...\n');
            [alignedData1, timeWindow1] = extract_aligned_traces(timeVector, Real_Signal_Smoothed, behavior1Times, preTime, postTime);
            if isempty(alignedData1)
                warning('No valid Behavior 1 trials after alignment in file %s.', csvFileName);
            else
                % Normalize the data to compute ΔF/F
                normalizedData1 = compute_deltaF_over_F(alignedData1, timeWindow1, preTime, postTime);

                % Append to aggregated data
                all_normalizedData1 = [all_normalizedData1; normalizedData1];
                if isempty(all_timeWindow1)
                    all_timeWindow1 = timeWindow1;
                end

                % Update global min and max for y-axis
                current_min1 = min(normalizedData1(:));
                current_max1 = max(normalizedData1(:));
                if current_min1 < global_min1
                    global_min1 = current_min1;
                end
                if current_max1 > global_max1
                    global_max1 = current_max1;
                end

                % Update global min and max for heatmap color scaling
                if current_min1 < global_heatmap_min1
                    global_heatmap_min1 = current_min1;
                end
                if current_max1 > global_heatmap_max1
                    global_heatmap_max1 = current_max1;
                end

                % Store normalized data in structure
                fileData(i).normalizedData1 = normalizedData1;
                fileData(i).timeWindow1 = timeWindow1;
            end
        else
            warning('No Behavior 1 tags found in Excel file %s.', xlsxFileName);
        end

        %% 10. Analyze Behavior2 Tags (First Pass)

        if ~isempty(behavior2Times)
            fprintf('Analyzing Behavior 2 tags...\n');
            [alignedData2, timeWindow2] = extract_aligned_traces(timeVector, Real_Signal_Smoothed, behavior2Times, preTime, postTime);
            if isempty(alignedData2)
                warning('No valid Behavior 2 trials after alignment in file %s.', csvFileName);
            else
                % Normalize the data to compute ΔF/F
                normalizedData2 = compute_deltaF_over_F(alignedData2, timeWindow2, preTime, postTime);

                % Append to aggregated data
                all_normalizedData2 = [all_normalizedData2; normalizedData2];
                if isempty(all_timeWindow2)
                    all_timeWindow2 = timeWindow2;
                end

                % Update global min and max for y-axis
                current_min2 = min(normalizedData2(:));
                current_max2 = max(normalizedData2(:));
                if current_min2 < global_min2
                    global_min2 = current_min2;
                end
                if current_max2 > global_max2
                    global_max2 = current_max2;
                end

                % Update global min and max for heatmap color scaling
                if current_min2 < global_heatmap_min2
                    global_heatmap_min2 = current_min2;
                end
                if current_max2 > global_heatmap_max2
                    global_heatmap_max2 = current_max2;
                end

                % Store normalized data in structure
                fileData(i).normalizedData2 = normalizedData2;
                fileData(i).timeWindow2 = timeWindow2;
            end
        else
            warning('No Behavior 2 tags found in Excel file %s.', xlsxFileName);
        end

        % Store file names
        fileData(i).csvFileName = csvFileName;
        fileData(i).xlsxFileName = xlsxFileName;
    end

    %% 11. Determine Global Y-Axis Limits with Padding

    padding = 0.1; % 10% padding

    % For Behavior1
    if isfinite(global_min1) && isfinite(global_max1)
        y_min1 = global_min1 - padding * abs(global_min1);
        y_max1 = global_max1 + padding * abs(global_max1);
    else
        y_min1 = -1; % Default values if no data
        y_max1 = 1;
        warning('No valid data found for Behavior 1. Using default y-axis limits.');
    end

    % For Behavior2
    if isfinite(global_min2) && isfinite(global_max2)
        y_min2 = global_min2 - padding * abs(global_min2);
        y_max2 = global_max2 + padding * abs(global_max2);
    else
        y_min2 = -1; % Default values if no data
        y_max2 = 1;
        warning('No valid data found for Behavior 2. Using default y-axis limits.');
    end

    %% 12. Second Pass: Generate Plots Using Collected Data and Global Limits

    for i = 1:numCsvFiles
        fprintf('\nGenerating Plots for File Pair %d of %d\n', i, numCsvFiles);

        % Retrieve stored data
        csvFileName = fileData(i).csvFileName;
        xlsxFileName = fileData(i).xlsxFileName;

        normalizedData1 = fileData(i).normalizedData1;
        timeWindow1 = fileData(i).timeWindow1;

        normalizedData2 = fileData(i).normalizedData2;
        timeWindow2 = fileData(i).timeWindow2;

        %% Generate Plots for Behavior1

        if ~isempty(normalizedData1)
            fprintf('Generating Plots for Behavior 1...\n');

            % Compute mean and SEM
            [meanTrace1, semTrace1] = compute_mean_sem(normalizedData1);

            % Figure for Behavior 1: Overlay Plot with Mean + SEM
            figure('Name', sprintf('Behavior 1: \DeltaF/F Signal Aligned to Onset - %s', csvFileName), 'NumberTitle', 'off');
            hold on;
            % Plot all individual trials
            plot(timeWindow1, normalizedData1', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
            % Plot mean trace
            plot(timeWindow1, meanTrace1, 'r', 'LineWidth', 2);
            % Plot SEM
            x_fill = [timeWindow1(:); flipud(timeWindow1(:))];
            y_fill = [meanTrace1(:) + semTrace1(:); flipud(meanTrace1(:) - semTrace1(:))];
            fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('\DeltaF/F');
            title(sprintf('Behavior 1: \DeltaF/F Signal Aligned to Onset - %s', csvFileName));
            legend({'Individual Trials', 'Mean', 'SEM'}, 'Location', 'best');
            grid on;
            hold off;
            % Set consistent y-axis limits
            ylim([y_min1, y_max1]);

            % Figure for Behavior 1 Heatmap
            figure('Name', sprintf('Behavior 1: Heatmap of \DeltaF/F Signal Aligned to Onset - %s', csvFileName), 'NumberTitle', 'off');
            imagesc(timeWindow1, 1:size(normalizedData1,1), normalizedData1);
            set(gca, 'YDir', 'normal');
            colormap('hot');
            colorbar;
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('Trial');
            title(sprintf('Behavior 1: Heatmap of \DeltaF/F Signal Aligned to Onset - %s', csvFileName));
            % Set consistent color scale
            caxis([global_heatmap_min1, global_heatmap_max1]);
        else
            warning('No normalized Behavior 1 data available for file %s.', csvFileName);
        end

        %% Generate Plots for Behavior2

        if ~isempty(normalizedData2)
            fprintf('Generating Plots for Behavior 2...\n');

            % Compute mean and SEM
            [meanTrace2, semTrace2] = compute_mean_sem(normalizedData2);

            % Figure for Behavior 2: Overlay Plot with Mean + SEM
            figure('Name', sprintf('Behavior 2: \DeltaF/F Signal Aligned to Onset - %s', csvFileName), 'NumberTitle', 'off');
            hold on;
            % Plot all individual trials
            plot(timeWindow2, normalizedData2', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
            % Plot mean trace
            plot(timeWindow2, meanTrace2, 'b', 'LineWidth', 2);
            % Plot SEM
            x_fill = [timeWindow2(:); flipud(timeWindow2(:))];
            y_fill = [meanTrace2(:) + semTrace2(:); flipud(meanTrace2(:) - semTrace2(:))];
            fill(x_fill, y_fill, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('\DeltaF/F');
            title(sprintf('Behavior 2: \DeltaF/F Signal Aligned to Onset - %s', csvFileName));
            legend({'Individual Trials', 'Mean', 'SEM'}, 'Location', 'best');
            grid on;
            hold off;
            % Set consistent y-axis limits
            ylim([y_min2, y_max2]);

            % Figure for Behavior 2 Heatmap
            figure('Name', sprintf('Behavior 2: Heatmap of \DeltaF/F Signal Aligned to Onset - %s', csvFileName), 'NumberTitle', 'off');
            imagesc(timeWindow2, 1:size(normalizedData2,1), normalizedData2);
            set(gca, 'YDir', 'normal');
            colormap('hot');
            colorbar;
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('Trial');
            title(sprintf('Behavior 2: Heatmap of \DeltaF/F Signal Aligned to Onset - %s', csvFileName));
            % Set consistent color scale
            caxis([global_heatmap_min2, global_heatmap_max2]);
        else
            warning('No normalized Behavior 2 data available for file %s.', csvFileName);
        end
    end

    %% 13. Generate Aggregate Figures for All Behavior1 and Behavior2 Data

    % Note: This section is outside the main loops to ensure it runs after all file pairs are processed.

    % Aggregate Figures for Behavior1
    if ~isempty(all_normalizedData1)
        fprintf('\nGenerating Aggregate Figures for Behavior 1...\n');

        % Define the AUC window: 0 to 100 seconds post tag
        aucStart = 0;
        aucEnd = 100;

        % Find indices corresponding to the AUC window
        auc_indices1 = all_timeWindow1 >= aucStart & all_timeWindow1 <= aucEnd;

        if sum(auc_indices1) == 0
            warning('No data points found in the AUC window (0 to 100 seconds) for Behavior 1. Skipping sorting.');
            sorted_normalizedData1 = all_normalizedData1;
        else
            % Compute AUC for each trial within the AUC window
            auc_values1 = trapz(all_timeWindow1(auc_indices1), all_normalizedData1(:, auc_indices1), 2);

            % Sort trials based on AUC in ascending order
            [sorted_auc1, sort_order1] = sort(auc_values1, 'ascend');

            % Rearrange the normalized data based on the sorted order
            sorted_normalizedData1 = all_normalizedData1(sort_order1, :);
        end

        % Compute overall mean and SEM for sorted data
        overall_mean1 = mean(sorted_normalizedData1, 1, 'omitnan');
        overall_sem1 = std(sorted_normalizedData1, 0, 1, 'omitnan') ./ sqrt(size(sorted_normalizedData1,1));

        % Aggregate Overlay Plot for Behavior1
        figure('Name', 'Aggregate Behavior 1: \DeltaF/F Signal Aligned to Onset (Sorted by AUC)', 'NumberTitle', 'off');
        hold on;
        % Plot all individual trials
        plot(all_timeWindow1, sorted_normalizedData1', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
        % Plot overall mean trace
        plot(all_timeWindow1, overall_mean1, 'k', 'LineWidth', 2);
        % Plot SEM
        x_fill = [all_timeWindow1(:); flipud(all_timeWindow1(:))];
        y_fill = [overall_mean1(:) + overall_sem1(:); flipud(overall_mean1(:) - overall_sem1(:))];
        fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xlabel('Time Relative to Behavior Onset (s)');
        ylabel('\DeltaF/F');
        title('Aggregate Behavior 1: \DeltaF/F Signal Aligned to Onset (Sorted by AUC)');
        legend({'Individual Trials', 'Mean', 'SEM'}, 'Location', 'best');
        grid on;
        hold off;
        % Set consistent y-axis limits
        ylim([y_min1, y_max1]);

        % Aggregate Sorted Heatmap for Behavior1
        figure('Name', 'Aggregate Behavior 1: Sorted Heatmap of \DeltaF/F Signal Aligned to Onset', 'NumberTitle', 'off');
        imagesc(all_timeWindow1, 1:size(sorted_normalizedData1,1), sorted_normalizedData1);
        set(gca, 'YDir', 'normal');
        colormap('hot');
        colorbar;
        xlabel('Time Relative to Behavior Onset (s)');
        ylabel('Trial (Sorted by AUC)');
        title('Aggregate Behavior 1: Sorted Heatmap of \DeltaF/F Signal Aligned to Onset');
        % Set consistent color scale
        caxis([global_heatmap_min1, global_heatmap_max1]);
    else
        warning('No aggregated Behavior 1 data available for aggregate figures.');
    end

    % Aggregate Figures for Behavior2
    if ~isempty(all_normalizedData2)
        fprintf('Generating Aggregate Figures for Behavior 2...\n');

        % Define the AUC window: 0 to 100 seconds post tag
        aucStart = 0;
        aucEnd = 100;

        % Find indices corresponding to the AUC window
        auc_indices2 = all_timeWindow2 >= aucStart & all_timeWindow2 <= aucEnd;

        if sum(auc_indices2) == 0
            warning('No data points found in the AUC window (0 to 100 seconds) for Behavior 2. Skipping sorting.');
            sorted_normalizedData2 = all_normalizedData2;
        else
            % Compute AUC for each trial within the AUC window
            auc_values2 = trapz(all_timeWindow2(auc_indices2), all_normalizedData2(:, auc_indices2), 2);

            % Sort trials based on AUC in ascending order
            [sorted_auc2, sort_order2] = sort(auc_values2, 'ascend');

            % Rearrange the normalized data based on the sorted order
            sorted_normalizedData2 = all_normalizedData2(sort_order2, :);
        end

        % Compute overall mean and SEM for sorted data
        overall_mean2 = mean(sorted_normalizedData2, 1, 'omitnan');
        overall_sem2 = std(sorted_normalizedData2, 0, 1, 'omitnan') ./ sqrt(size(sorted_normalizedData2,1));

        % Aggregate Overlay Plot for Behavior2
        figure('Name', 'Aggregate Behavior 2: \DeltaF/F Signal Aligned to Onset (Sorted by AUC)', 'NumberTitle', 'off');
        hold on;
        % Plot all individual trials
        plot(all_timeWindow2, sorted_normalizedData2', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
        % Plot overall mean trace
        plot(all_timeWindow2, overall_mean2, 'k', 'LineWidth', 2);
        % Plot SEM
        x_fill = [all_timeWindow2(:); flipud(all_timeWindow2(:))];
        y_fill = [overall_mean2(:) + overall_sem2(:); flipud(overall_mean2(:) - overall_sem2(:))];
        fill(x_fill, y_fill, 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        xlabel('Time Relative to Behavior Onset (s)');
        ylabel('\DeltaF/F');
        title('Aggregate Behavior 2: \DeltaF/F Signal Aligned to Onset (Sorted by AUC)');
        legend({'Individual Trials', 'Mean', 'SEM'}, 'Location', 'best');
        grid on;
        hold off;
        % Set consistent y-axis limits
        ylim([y_min2, y_max2]);

        % Aggregate Sorted Heatmap for Behavior2
        figure('Name', 'Aggregate Behavior 2: Sorted Heatmap of \DeltaF/F Signal Aligned to Onset', 'NumberTitle', 'off');
        imagesc(all_timeWindow2, 1:size(sorted_normalizedData2,1), sorted_normalizedData2);
        set(gca, 'YDir', 'normal');
        colormap('hot');
        colorbar;
        xlabel('Time Relative to Behavior Onset (s)');
        ylabel('Trial (Sorted by AUC)');
        title('Aggregate Behavior 2: Sorted Heatmap of \DeltaF/F Signal Aligned to Onset');
        % Set consistent color scale
        caxis([global_heatmap_min2, global_heatmap_max2]);
    else
        warning('No aggregated Behavior 2 data available for aggregate figures.');
    end

    %% 14. Helper Functions

    function [alignedData, timeWindow] = extract_aligned_traces(timeVec, signal, tagTimes, pre, post)
        % Extracts signal segments aligned to each tag time.
        % Inputs:
        %   timeVec - vector of time points in seconds
        %   signal - vector of Real_Signal_Smoothed values
        %   tagTimes - vector of behavior onset times in seconds
        %   pre - time before tag to include (negative value)
        %   post - time after tag to include
        % Outputs:
        %   alignedData - matrix of size (numTags x numTimePoints)
        %   timeWindow - vector of relative time points

        numTags = length(tagTimes);
        samplingInterval = median(diff(timeVec));  % Estimate the sampling interval
        numTimePoints = round((post - pre) / samplingInterval) + 1;
        timeWindow = linspace(pre, post, numTimePoints)';

        % Initialize aligned data matrix with NaNs
        alignedData = NaN(numTags, length(timeWindow));

        for i = 1:numTags
            tagTime = tagTimes(i);
            % Define desired absolute time points for the window
            desiredTimes = tagTime + timeWindow;

            % Check if desiredTimes are within the range of timeVec
            if desiredTimes(1) < timeVec(1) || desiredTimes(end) > timeVec(end)
                warning('Trial %d: Window exceeds data limits. Trial will be excluded.', i);
                continue;  % Skip this trial
            end

            % Interpolate signal at desired times
            interpolatedSignal = interp1(timeVec, signal, desiredTimes, 'linear', NaN);

            % Assign to aligned data if no NaNs
            if all(isfinite(interpolatedSignal))
                alignedData(i, :) = interpolatedSignal';
            else
                warning('Trial %d: Interpolated signal contains NaNs. Trial will be excluded.', i);
            end
        end

        % Remove trials with any NaNs
        validTrials = all(isfinite(alignedData), 2);
        if any(~validTrials)
            numExcluded = sum(~validTrials);
            warning('%d trials have incomplete data and will be excluded.', numExcluded);
            alignedData = alignedData(validTrials, :);
        end
    end

    function [meanTrace, semTrace] = compute_mean_sem(alignedData)
        % Computes mean and SEM across trials.
        % Inputs:
        %   alignedData - matrix of size (numTrials x numTimePoints)
        % Outputs:
        %   meanTrace - mean across trials
        %   semTrace - standard error of the mean across trials

        meanTrace = mean(alignedData, 1, 'omitnan');
        semTrace = std(alignedData, 0, 1, 'omitnan') ./ sqrt(size(alignedData,1));
    end

    function normalizedData = compute_deltaF_over_F(alignedData, timeWindow, preTime, postTime)
        % Computes Delta F/F for aligned data.
        % Inputs:
        %   alignedData - matrix of size (numTrials x numTimePoints)
        %   timeWindow  - vector of relative time points
        %   preTime     - start time of the window (negative value)
        %   postTime    - end time of the window (positive value)
        % Output:
        %   normalizedData - matrix of Delta F/F values

        % Identify indices corresponding to the baseline period (-5 to 0 seconds)
        baselineStart = -5; % seconds before tag
        baselineEnd = 0;    % seconds relative to tag
        baselineIdx = timeWindow >= baselineStart & timeWindow <= baselineEnd;

        if sum(baselineIdx) == 0
            error('No data points found in the baseline period (-5 to 0 seconds).');
        end

        % Compute F_baseline for each trial
        F_baseline = mean(alignedData(:, baselineIdx), 2);  % (numTrials x 1)

        % Handle trials with zero or near-zero baseline to avoid division by zero
        zeroBaseline = F_baseline < 1e-6;
        if any(zeroBaseline)
            warning('%d trials have a baseline F < 1e-6 and will be excluded from normalization.', sum(zeroBaseline));
            F_baseline(zeroBaseline) = NaN;  % Assign NaN to exclude these trials
        end

        % Compute Delta F/F
        normalizedData = (alignedData - F_baseline) ./ F_baseline;  % (numTrials x numTimePoints)

        % Exclude trials with NaN baseline
        validTrials = all(isfinite(normalizedData), 2);
        if any(~validTrials)
            numExcluded = sum(~validTrials);
            warning('%d trials have invalid normalization and will be excluded.', numExcluded);
            normalizedData = normalizedData(validTrials, :);
        end
    end

end
