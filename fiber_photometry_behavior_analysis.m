% fiber_photometry_behavior_analysis.m
% This script analyzes and visualizes behavior-tagged fiber photometry data.
% It performs the following:
%   1. Imports raw CSV data and behavior tags from an Excel file.
%   2. Applies a low-pass Butterworth filter to smooth the Real_Signal.
%   3. Aligns the smoothed signal to behavior onset times.
%   4. Computes ΔF/F normalization.
%   5. Generates overlay plots with mean and SEM.
%   6. Creates heatmaps of the normalized signals.

function fiber_photometry_behavior_analysis()
    % Clear workspace and close all figures
    close all;
    clearvars;
    
    %% 1. Import Data
    
    % Prompt user to select the CSV raw data file
    [csvFileName, csvPath] = uigetfile('*.csv', 'Select the Raw CSV Data File');
    if isequal(csvFileName,0)
        disp('User canceled CSV file selection.');
        return;
    end
    csvFullPath = fullfile(csvPath, csvFileName);
    fprintf('Loading CSV data from: %s\n', csvFullPath);
    
    try
        rawData = readtable(csvFullPath);
    catch ME
        error('Failed to read CSV file: %s', ME.message);
    end
    
    % Check if required columns exist
    requiredColumns = [1, 2, 4];
    if width(rawData) < max(requiredColumns)
        error('CSV file does not contain the required columns (at least 4 columns expected).');
    end
    
    % Extract timestamps and LED signals
    try
        TimeStamp = rawData{:,1};  % Column 1
        LED_410 = rawData{:,2};    % Column 2
        LED_470 = rawData{:,4};    % Column 4
    catch ME
        error('Error extracting data from CSV: %s', ME.message);
    end

    % Validate and convert TimeStamp to seconds
    if isduration(TimeStamp)
        timeVector = seconds(TimeStamp) - seconds(TimeStamp(1));
    elseif isdatetime(TimeStamp)
        timeVector = seconds(TimeStamp - TimeStamp(1));
    else
        error('Timestamp column must be of type datetime or duration.');
    end
    timeVector = timeVector(:);  % Ensure column vector

    % Calculate Real_Signal (470/410)
    Real_Signal = LED_470 ./ LED_410;
    Real_Signal = Real_Signal(:);  % Ensure column vector

    % Handle non-finite values
    finiteIdx = isfinite(timeVector) & isfinite(Real_Signal);
    if any(~finiteIdx)
        warning('Non-finite values found in the data. These will be removed.');
        timeVector = timeVector(finiteIdx);
        Real_Signal = Real_Signal(finiteIdx);
    end

    % Check if data is sufficient after cleaning
    if isempty(timeVector) || isempty(Real_Signal)
        error('No valid data available after cleaning.');
    end

    %% 2. Apply Low-Pass Butterworth Filter to Smooth the Signal
    
    % Define filter parameters
    % These parameters can be adjusted based on your data's characteristics
    % For example, Fc = 5 Hz is a common choice, but you may need to adjust it
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
    
    % Optional: Plot raw vs. smoothed signal for verification
    figure('Name', 'Raw vs. Smoothed Real_Signal', 'NumberTitle', 'off');
    plot(timeVector, Real_Signal, 'b', 'DisplayName', 'Raw Signal');
    hold on;
    plot(timeVector, Real_Signal_Smoothed, 'r', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Signal');
    xlabel('Time (s)');
    ylabel('470/410 Signal');
    title('Raw vs. Smoothed Real\_Signal');
    legend('show');
    grid on;
    hold off;
    
    %% 3. Import Behavior Tags
    
    % Prompt user to select the BehaviorTags Excel file
    [excelFileName, excelPath] = uigetfile('*.xlsx', 'Select the Behavior Tags Excel File');
    if isequal(excelFileName,0)
        disp('User canceled Excel file selection.');
        return;
    end
    excelFullPath = fullfile(excelPath, excelFileName);
    fprintf('Loading Behavior Tags from: %s\n', excelFullPath);
    
    % Read sheets from Excel file
    try
        % Assuming sheets are named 'Behavior1' and 'Behavior2'
        behavior1Data = readtable(excelFullPath, 'Sheet', 'Behavior1');
        behavior2Data = readtable(excelFullPath, 'Sheet', 'Behavior2');
    catch ME
        error('Failed to read Behavior Tags Excel file: %s', ME.message);
    end
    
    % Extract behavior onset times
    if ismember('Behavior1_Onset_Time_s', behavior1Data.Properties.VariableNames)
        behavior1Times = behavior1Data.Behavior1_Onset_Time_s;
    else
        error('Sheet "Behavior1" does not contain the column "Behavior1_Onset_Time_s".');
    end
    
    if ismember('Behavior2_Onset_Time_s', behavior2Data.Properties.VariableNames)
        behavior2Times = behavior2Data.Behavior2_Onset_Time_s;
    else
        error('Sheet "Behavior2" does not contain the column "Behavior2_Onset_Time_s".');
    end

    %% 4. Define Analysis Parameters
    
    % Define time window around behavior onset
    preTime = -30;  % seconds before tag
    postTime = 30; % seconds after tag
    windowDuration = postTime - preTime; % total window duration
    
    %% 5. Analyze and Visualize Behavior1 Tags
    
    if ~isempty(behavior1Times)
        fprintf('Analyzing Behavior 1 tags...\n');
        [alignedData1, timeWindow] = extract_aligned_traces(timeVector, Real_Signal_Smoothed, behavior1Times, preTime, postTime);
        if isempty(alignedData1)
            warning('No valid Behavior 1 trials after alignment.');
        else
            % Normalize the data to compute ΔF/F
            normalizedData1 = compute_deltaF_over_F(alignedData1, timeWindow, preTime, postTime);
            
            [meanTrace1, semTrace1] = compute_mean_sem(normalizedData1);
            
            % Figure 1 for Behavior 1: Overlay Plot with Mean + SEM
            figure('Name', 'Behavior 1: ΔF/F Signal Aligned to Onset', 'NumberTitle', 'off');
            hold on;
            % Plot all individual trials
            plot(timeWindow, normalizedData1', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
            % Plot mean trace
            plot(timeWindow, meanTrace1, 'r', 'LineWidth', 2);
            % Plot SEM
            % Ensure both X and Y are column vectors
            x_fill = [timeWindow(:); flipud(timeWindow(:))];
            y_fill = [meanTrace1(:) + semTrace1(:); flipud(meanTrace1(:) - semTrace1(:))];
            fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('\DeltaF/F');
            title('Behavior 1: \DeltaF/F Signal Aligned to Onset');
            legend({'Individual Trials', 'Mean', 'SEM'}, 'Location', 'best');
            grid on;
            hold off;
            
            % Figure 2 for Behavior 1 Heatmap
            figure('Name', 'Behavior 1: Heatmap of ΔF/F Signal Aligned to Onset', 'NumberTitle', 'off');
            imagesc(timeWindow, 1:size(normalizedData1,1), normalizedData1);
            set(gca, 'YDir', 'normal');
            colormap('jet');
            colorbar;
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('Trial');
            title('Behavior 1: Heatmap of \DeltaF/F Signal Aligned to Onset');
        end
    else
        warning('No Behavior 1 tags found.');
    end

    %% 6. Analyze and Visualize Behavior2 Tags
    
    if ~isempty(behavior2Times)
        fprintf('Analyzing Behavior 2 tags...\n');
        [alignedData2, timeWindow] = extract_aligned_traces(timeVector, Real_Signal_Smoothed, behavior2Times, preTime, postTime);
        if isempty(alignedData2)
            warning('No valid Behavior 2 trials after alignment.');
        else
            % Normalize the data to compute ΔF/F
            normalizedData2 = compute_deltaF_over_F(alignedData2, timeWindow, preTime, postTime);
            
            [meanTrace2, semTrace2] = compute_mean_sem(normalizedData2);
            
            % Figure 3 for Behavior 2: Overlay Plot with Mean + SEM
            figure('Name', 'Behavior 2: ΔF/F Signal Aligned to Onset', 'NumberTitle', 'off');
            hold on;
            % Plot all individual trials
            plot(timeWindow, normalizedData2', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
            % Plot mean trace
            plot(timeWindow, meanTrace2, 'b', 'LineWidth', 2);
            % Plot SEM
            % Ensure both X and Y are column vectors
            x_fill = [timeWindow(:); flipud(timeWindow(:))];
            y_fill = [meanTrace2(:) + semTrace2(:); flipud(meanTrace2(:) - semTrace2(:))];
            fill(x_fill, y_fill, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('\DeltaF/F');
            title('Behavior 2: \DeltaF/F Signal Aligned to Onset');
            legend({'Individual Trials', 'Mean', 'SEM'}, 'Location', 'best');
            grid on;
            hold off;
            
            % Figure 4 for Behavior 2 Heatmap
            figure('Name', 'Behavior 2: Heatmap of ΔF/F Signal Aligned to Onset', 'NumberTitle', 'off');
            imagesc(timeWindow, 1:size(normalizedData2,1), normalizedData2);
            set(gca, 'YDir', 'normal');
            colormap('jet');
            colorbar;
            xlabel('Time Relative to Behavior Onset (s)');
            ylabel('Trial');
            title('Behavior 2: Heatmap of \DeltaF/F Signal Aligned to Onset');
        end
    else
        warning('No Behavior 2 tags found.');
    end

    %% 7. Helper Functions
    
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
        baselineIdx = timeWindow >= preTime & timeWindow <= 0;

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
