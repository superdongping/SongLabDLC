function fiber_photometry_video_sync_v4_3  
clc;
warning off;
close all;    
    % Create a UI figure
    fig = uifigure('Name', 'Fiber Photometry and Video Synchronization', ...
                   'Position', [50, 50, 1600, 1000]); % Increased width to accommodate new UI elements
    
    %% Existing UI Components
    
    % Add a button to load the CSV file
    csvButton = uibutton(fig, 'Text', 'Load CSV', ...
        'Position', [20, 50, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) loadCSV());
    
    % Add a button to load the AVI file (initially disabled)
    aviButton = uibutton(fig, 'Text', 'Load Video', ...
        'Position', [140, 50, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) loadAVI());
    aviButton.Enable = 'off'; % Enable only after CSV is loaded successfully
    
    % Add "Play" and "Pause" buttons
    playButton = uibutton(fig, 'Text', 'Play', ...
        'Position', [260, 50, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) playVideo());
    pauseButton = uibutton(fig, 'Text', 'Pause', ...
        'Position', [380, 50, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) pauseVideo());
    
    % Add "Tag Behavior 1" and "Tag Behavior 2" buttons
    tagBehavior1Button = uibutton(fig, 'Text', 'Tag Behavior 1', ...
        'Position', [500, 50, 120, 30], ...
        'ButtonPushedFcn', @(btn, event) tagBehavior(1));
    
    tagBehavior2Button = uibutton(fig, 'Text', 'Tag Behavior 2', ...
        'Position', [640, 50, 120, 30], ...
        'ButtonPushedFcn', @(btn, event) tagBehavior(2));
    
    % Add "Save Tags" button
    saveTagsButton = uibutton(fig, 'Text', 'Save Tags', ...
        'Position', [780, 50, 100, 30], ...
        'ButtonPushedFcn', @(btn, event) saveTags(), ...
        'Enable', 'off'); % Disabled until tags are created
    
    %% *** New Buttons for Time Navigation ***
    
    % Add "5 Second Backward" button
    backwardButton = uibutton(fig, 'Text', '5 sec Backward', ...
        'Position', [900, 50, 120, 30], ...
        'ButtonPushedFcn', @(btn, event) moveTimeBackward());
    
    % Add "5 Second Forward" button
    forwardButton = uibutton(fig, 'Text', '5 sec Forward', ...
        'Position', [1030, 50, 120, 30], ...
        'ButtonPushedFcn', @(btn, event) moveTimeForward());
    
    %% *** End of New Buttons ***
    
    %% *** New UI Components ***
    
    % Add "Baseline Correction" checkbox
    baselineCorrectionCheckbox = uicheckbox(fig, 'Text', 'Baseline Correction', ...
        'Position', [1160, 55, 150, 20], ...
        'ValueChangedFcn', @(cb, event) toggleBaselineCorrection(cb.Value));
    
    % Add "Subtract Low-Pass Filter" checkbox
    subtractLowPassCheckbox = uicheckbox(fig, 'Text', 'Subtract Low-Pass Filter', ...
        'Position', [1160, 25, 180, 20], ... % Positioned below the Baseline Correction checkbox
        'ValueChangedFcn', @(cb, event) toggleSubtractLowPass(cb.Value));
    
    % Add "Cutoff Frequency (Hz)" label
    cutoffFreqLabel = uilabel(fig, 'Text', 'Cutoff Frequency (Hz):', ...
        'Position', [1160, 85, 150, 20], ...
        'HorizontalAlignment', 'left');
    
    % Add "Cutoff Frequency (Hz)" input box
    cutoffFreqInput = uieditfield(fig, 'numeric', ...
        'Position', [1310, 80, 100, 30], ...
        'Value', 0.05, ... % Default cutoff frequency
        'Limits', [0.001, 100], ... % Reasonable range for cutoff frequency
        'ValueChangedFcn', @(ef, event) cutoffFreqChanged(ef.Value));
    
    %% *** End of New UI Components ***
    
    % Add a slider to scroll through time
    timeSlider = uislider(fig, 'Position', [60, 150, 1350, 20], ... % Adjusted width to match increased figure size
        'ValueChangedFcn', @(sld, event) updateTime());
    timeSlider.Limits = [0, 1];  % Initialize with 0-1, will be updated after loading data
    
    % Display current time label
    timeLabel = uilabel(fig, 'Position', [700, 120, 200, 20], 'Text', 'Time: 0.00 s', ...
        'FontSize', 12, 'HorizontalAlignment', 'center'); % Adjusted position for larger figure
    
    % Axes for displaying the video
    axVideo = uiaxes(fig, 'Position', [100, 420, 900, 600]);
    axis(axVideo, 'off');  % Hide axes for better video display
    
    % Axes for plotting 470/410 signal
    axSignal = uiaxes(fig, 'Position', [20, 180, 1400, 250]);
    title(axSignal, '470/410 Signal');
    xlabel(axSignal, 'Time (s)');
    ylabel(axSignal, '470/410');
    hold(axSignal, 'on'); % Hold on to add dynamic elements
    
    %% Variables to store data
    csvData = [];
    timeVector = [];
    Real_Signal = [];
    Real_Signal_Smoothed = [];
    aviObj = [];
    videoFrameRate = 30;  % Default frame rate, will be updated after loading video
    videoTimer = timer('ExecutionMode', 'fixedRate', ...
        'Period', 1/videoFrameRate, ...
        'TimerFcn', @(~,~) playFrame());
    interpolatedTime = [];
    interpolatedSignal = [];
    currentTimeLine = []; % Handle for the current time indicator line
    
    %% Variables to store behavior tags
    behavior1Times = []; % Stores times for Behavior 1
    behavior2Times = []; % Stores times for Behavior 2
    csvBaseName = '';    % Stores the base name of the loaded CSV file
    csvPath = '';        % Stores the path of the loaded CSV file
    
    %% *** New Variables ***
    baselineCorrectionEnabled = false;        % Existing variable for baseline correction
    subtractLowPassEnabled = false;          % New variable for low-pass subtraction
    cutoffFrequency = 0.05;                   % Default cutoff frequency (Hz)
    %% *** End of New Variables ***
    
    %% Callback Functions
    
    function loadCSV()
        % Let the user select a CSV file
        [csvFileName, csvPathSelected] = uigetfile('*.csv', 'Select the CSV file');
        if isequal(csvFileName,0)
            return; % If the user cancels, do nothing
        end
        csvPath = csvPathSelected; % Store the path for saving tags
        try
            csvData = readtable(fullfile(csvPath, csvFileName));
        catch ME
            uialert(fig, ['Failed to read CSV file: ' ME.message], 'Error');
            return;
        end
        
        % Store the base name without extension
        [~, csvBaseName, ~] = fileparts(csvFileName);
        
        % Check if required columns exist
        % Updated to check for at least 3 or 4 columns depending on the data
        if width(csvData) < 3
            uialert(fig, 'CSV file does not contain the required columns (at least 3 columns expected).', 'Error');
            return;
        end
        
        % Extract timestamps and LED signals
        try
            TimeStamp = csvData{:, 1};  % Extracted timestamps as datetime or duration type
            LED_410 = csvData{:, 2};
            
            % *** Modified Section Starts Here ***
            % Dynamically identify the LED_470 column (either column 3 or 4)
            colNames = csvData.Properties.VariableNames;
            index470 = find(contains(lower(colNames), '470'));
            
            if ~isempty(index470)
                % If a column name contains '470', use it
                LED_470 = csvData{:, index470(1)};
            elseif width(csvData) >=4
                % If no '470' in column names but at least 4 columns, use column 4
                LED_470 = csvData{:, 4};
            else
                % Fallback: use column 3
                LED_470 = csvData{:, 3};
                uialert(fig, '470 data assumed to be in column 3 as no specific identifier was found.', 'Warning');
            end
            % *** Modified Section Ends Here ***
            
        catch ME
            uialert(fig, ['Error extracting data from CSV: ' ME.message], 'Error');
            return;
        end

        % Validate TimeStamp
        if ~isdatetime(TimeStamp) && ~isduration(TimeStamp)
            uialert(fig, 'Timestamp column must be datetime or duration type.', 'Error');
            return;
        end

        % Convert duration or datetime to seconds and ensure it's a column vector
        if isduration(TimeStamp)
            timeVector = seconds(TimeStamp) - seconds(TimeStamp(1));
        else
            timeVector = seconds(TimeStamp - TimeStamp(1));
        end
        timeVector = timeVector(:);  % Ensure it's a column vector

        % Calculate the real signal (470/410) and ensure it's a column vector
        Real_Signal = LED_470 ./ LED_410;
        Real_Signal = Real_Signal(:);  % Ensure it's a column vector

        % Handle non-finite values in Real_Signal
        finiteIdx = isfinite(timeVector) & isfinite(Real_Signal);
        if any(~finiteIdx)
            warning('Non-finite values found in the data. These will be removed.');
            timeVector = timeVector(finiteIdx);
            Real_Signal = Real_Signal(finiteIdx);
        end

        % Check if data is sufficient after cleaning
        if isempty(timeVector) || isempty(Real_Signal)
            uialert(fig, 'No valid data available after cleaning.', 'Error');
            return;
        end

        % Read the current cutoff frequency from the input box
        cutoffFrequency = cutoffFreqInput.Value;
        if isempty(cutoffFrequency) || ~isnumeric(cutoffFrequency) || cutoffFrequency <= 0
            uialert(fig, 'Please enter a valid positive number for the cutoff frequency.', 'Invalid Input');
            cutoffFreqInput.Value = 0.05; % Reset to default
            cutoffFrequency = 0.05;
        end

        % Apply Low-Pass Butterworth Filter to Smooth the Signal
        filterOrder = 4;    % Order of the Butterworth filter

        % Estimate the sampling frequency (Fs)
        samplingInterval = median(diff(timeVector));  % Median sampling interval in seconds
        Fs = 1 / samplingInterval;                     % Sampling frequency in Hz
        fprintf('Estimated Sampling Frequency: %.2f Hz\n', Fs);

        % Validate that cutoffFrequency is less than Nyquist frequency
        if cutoffFrequency >= Fs/2
            uialert(fig, ['Cutoff frequency must be less than half the sampling rate (Fs/2 = ', num2str(Fs/2), ' Hz).'], 'Invalid Cutoff Frequency');
            cutoffFreqInput.Value = min(cutoffFrequency, Fs/2 - 0.01); % Adjust to just below Nyquist
            cutoffFrequency = cutoffFreqInput.Value;
            fprintf('Adjusted Cutoff Frequency to %.4f Hz\n', cutoffFrequency);
        end

        % Design Butterworth low-pass filter
        [b, a] = butter(filterOrder, cutoffFrequency/(Fs/2), 'low');

        % Apply zero-phase filtering using filtfilt to prevent phase distortion
        Real_Signal_Smoothed = filtfilt(b, a, Real_Signal);
        
        %% *** Apply Baseline Correction if Enabled ***
        if baselineCorrectionEnabled
            Real_Signal_Smoothed = applyBaselineCorrection(Real_Signal_Smoothed, timeVector);
        end
        %% *** End of Baseline Correction ***
        
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

        % Plot the smoothed (and possibly corrected) signal on axSignal
        processAndPlotSignal();
        
        % Update the slider range based on the time vector
        timeSlider.Limits = [0, max(timeVector)];
        timeSlider.Value = 0;  % Reset slider to start

        % Update the time label
        timeLabel.Text = sprintf('Time: %.2f s', timeSlider.Value);

        % Enable the Load Video button now that CSV is loaded
        aviButton.Enable = 'on';

        % Initialize the current time indicator line
        if isgraphics(currentTimeLine)
            delete(currentTimeLine);
        end
        currentTimeLine = plot(axSignal, [0 0], ylim(axSignal), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Current Time');

        % Reset behavior tags when a new CSV is loaded
        behavior1Times = [];
        behavior2Times = [];
        saveTagsButton.Enable = 'off';
    end

    function loadAVI()
        % Ensure that CSV data is loaded
        if isempty(timeVector) || isempty(Real_Signal_Smoothed)
            uialert(fig, 'Please load a valid CSV file first.', 'Error');
            return;
        end

        % Let the user select an AVI file
        [aviFileName, aviPathSelected] = uigetfile('*.avi', 'Select the AVI file');
        if isequal(aviFileName, 0)
            return; % If the user cancels, do nothing
        end
        aviPath = aviPathSelected; % Store path if needed

        try
            aviObj = VideoReader(fullfile(aviPath, aviFileName));
        catch ME
            uialert(fig, ['Failed to read AVI file: ' ME.message], 'Error');
            return;
        end
        videoFrameRate = aviObj.FrameRate;
        videoTimer.Period = 1 / videoFrameRate;  % Update timer period based on video frame rate

        % Interpolate the signal to match the video frame rate
        interpolatedTime = (0:1/videoFrameRate:max(timeVector))';
        
        %% *** Apply Baseline Correction to Interpolated Signal if Enabled ***
        if baselineCorrectionEnabled
            interpolatedSignal = applyBaselineCorrection(interp1(timeVector, Real_Signal_Smoothed, interpolatedTime, 'linear', 'extrap'), interpolatedTime);
        else
            interpolatedSignal = interp1(timeVector, Real_Signal_Smoothed, interpolatedTime, 'linear', 'extrap');
        end
        %% *** End of Baseline Correction ***
        
        % Check for finite values in interpolatedSignal
        if any(~isfinite(interpolatedSignal))
            uialert(fig, 'Interpolated signal contains non-finite values.', 'Error');
            return;
        end

        % Update the signal plot with interpolated data
        processAndPlotSignal();

        % Update the slider range based on the interpolated time
        timeSlider.Limits = [0, max(interpolatedTime)];
        timeSlider.Value = 0;  % Reset slider to start

        % Update the time label
        timeLabel.Text = sprintf('Time: %.2f s', timeSlider.Value);

        % Display the first frame of the video
        aviObj.CurrentTime = 0;
        try
            frame = readFrame(aviObj);
            imshow(frame, 'Parent', axVideo);
        catch ME
            uialert(fig, ['Failed to read the first frame of the video: ' ME.message], 'Error');
            return;
        end
    end

    function playVideo()
        % Start the timer to play the video from the current time
        if isempty(aviObj)
            uialert(fig, 'Please load a video first.', 'Error');
            return;
        end
        if strcmp(videoTimer.Running, 'off')
            start(videoTimer);
        end
    end

    function pauseVideo()
        % Stop the timer to pause the video
        stop(videoTimer);
    end

    function playFrame()
        if isempty(aviObj) || isempty(interpolatedTime)
            return; % Do nothing if data or video is not loaded
        end

        % Display the current frame of the video
        if hasFrame(aviObj)
            try
                frame = readFrame(aviObj);
                imshow(frame, 'Parent', axVideo);
            catch
                stop(videoTimer); % Stop if unable to read frame
                return;
            end

            % Update the slider and plot based on the current time
            currentTime = aviObj.CurrentTime;
            timeSlider.Value = currentTime;
            updatePlot(currentTime);
        else
            stop(videoTimer); % Stop the video when it reaches the end
        end
    end

    function updateTime()
        % Get the current time from the slider
        currentTime = timeSlider.Value;
        
        % Update the time label
        timeLabel.Text = sprintf('Time: %.2f s', currentTime);
        
        % Adjust the video to the selected time and update frame
        updateVideoFrame(currentTime);
        updatePlot(currentTime);
    end

    function updatePlot(currentTime)
        % Update the current time indicator line
        if isgraphics(currentTimeLine)
            currentTimeLine.XData = [currentTime currentTime];
        else
            currentTimeLine = plot(axSignal, [currentTime currentTime], ylim(axSignal), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Current Time');
        end
    end

    function updateVideoFrame(currentTime)
        % Update the video frame based on the current time
        if isempty(aviObj)
            return;
        end

        aviObj.CurrentTime = currentTime;
        try
            frame = readFrame(aviObj);
            imshow(frame, 'Parent', axVideo);
        catch
            % If unable to read frame, do nothing or handle accordingly
        end
    end

    %% *** New Callback Functions for Time Navigation ***
    
    function moveTimeForward()
        % Move the time slider forward by 5 seconds
        currentTime = timeSlider.Value;
        newTime = currentTime + 5;
        if newTime > timeSlider.Limits(2)
            newTime = timeSlider.Limits(2);
        end
        timeSlider.Value = newTime;
        updateTime();
    end

    function moveTimeBackward()
        % Move the time slider backward by 5 seconds
        currentTime = timeSlider.Value;
        newTime = currentTime - 5;
        if newTime < timeSlider.Limits(1)
            newTime = timeSlider.Limits(1);
        end
        timeSlider.Value = newTime;
        updateTime();
    end
    
    %% *** End of New Callback Functions ***
    
    %% *** New Callback Function: Baseline Correction Toggle ***
    function toggleBaselineCorrection(enable)
        baselineCorrectionEnabled = enable;
        
        % If data is already loaded, reprocess and update the plot
        if ~isempty(timeVector) && ~isempty(Real_Signal)
            processAndPlotSignal();
        end
    end
    %% *** End of New Callback Function ***
    
    %% *** New Callback Function: Subtract Low-Pass Filter Toggle ***
    function toggleSubtractLowPass(enable)
        subtractLowPassEnabled = enable;
        
        % If data is already loaded, reprocess and update the plot
        if ~isempty(timeVector) && ~isempty(Real_Signal)
            processAndPlotSignal();
        end
    end
    %% *** End of New Callback Function ***
    
    %% *** New Callback Function: Cutoff Frequency Changed ***
    function cutoffFreqChanged(newFc)
        % Validate the new cutoff frequency
        if isempty(newFc) || ~isnumeric(newFc) || newFc <= 0
            uialert(fig, 'Please enter a valid positive number for the cutoff frequency.', 'Invalid Input');
            cutoffFreqInput.Value = cutoffFrequency; % Revert to previous valid value
            return;
        end
        
        % Update the cutoff frequency
        cutoffFrequency = newFc;
        
        % Reprocess the signal with the new cutoff frequency
        if ~isempty(timeVector) && ~isempty(Real_Signal)
            processAndPlotSignal();
        end
    end
    %% *** End of New Callback Function ***
    
    %% *** New Function: Process and Plot Signal ***
    function processAndPlotSignal()
        % Start with the original signal
        processedSignal = Real_Signal;
        
        % Read the current cutoff frequency from the input box
        cutoffFrequency = cutoffFreqInput.Value;
        if isempty(cutoffFrequency) || ~isnumeric(cutoffFrequency) || cutoffFrequency <= 0
            uialert(fig, 'Please enter a valid positive number for the cutoff frequency.', 'Invalid Input');
            cutoffFreqInput.Value = 0.05; % Reset to default
            cutoffFrequency = 0.05;
        end
        
        % Estimate the sampling frequency (Fs)
        samplingInterval = median(diff(timeVector));  % Median sampling interval in seconds
        Fs = 1 / samplingInterval;                     % Sampling frequency in Hz
        fprintf('Estimated Sampling Frequency: %.2f Hz\n', Fs);
        
        % Validate that cutoffFrequency is less than Nyquist frequency
        if cutoffFrequency >= Fs/2
            uialert(fig, ['Cutoff frequency must be less than half the sampling rate (Fs/2 = ', num2str(Fs/2), ' Hz).'], 'Invalid Cutoff Frequency');
            cutoffFreqInput.Value = min(cutoffFrequency, Fs/2 - 0.01); % Adjust to just below Nyquist
            cutoffFrequency = cutoffFreqInput.Value;
            fprintf('Adjusted Cutoff Frequency to %.4f Hz\n', cutoffFrequency);
        end
        
        % Apply Low-Pass Butterworth Filter to Smooth the Signal
        filterOrder = 4;    % Order of the Butterworth filter
        
        % Design Butterworth low-pass filter
        [b, a] = butter(filterOrder, cutoffFrequency/(Fs/2), 'low');
        
        % Apply zero-phase filtering using filtfilt to prevent phase distortion
        Real_Signal_Smoothed = filtfilt(b, a, Real_Signal);
        
        % Apply baseline correction if enabled
        if baselineCorrectionEnabled
            Real_Signal_Smoothed = applyBaselineCorrection(Real_Signal_Smoothed, timeVector);
        end
        
        % Determine the processed signal based on user selections
        if subtractLowPassEnabled
            processedSignal = Real_Signal - Real_Signal_Smoothed;
        else
            processedSignal = Real_Signal_Smoothed;
        end
        
        % Update the plot on axSignal
        cla(axSignal); % Clear existing plot
        
        if subtractLowPassEnabled
            plot(axSignal, timeVector, processedSignal, 'g', 'LineWidth', 1.5, 'DisplayName', 'High-Frequency Signal');
        else
            plot(axSignal, timeVector, processedSignal, 'r', 'LineWidth', 1.5, 'DisplayName', 'Processed Signal');
        end
        
        xlim(axSignal, [0, max(timeVector)]);
        
        % Update the current time indicator line
        if isgraphics(currentTimeLine)
            currentTimeLine.XData = [timeSlider.Value timeSlider.Value];
        else
            currentTimeLine = plot(axSignal, [timeSlider.Value timeSlider.Value], ylim(axSignal), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Current Time');
        end
        
        % If AVI is loaded, re-interpolate with the new signal
        if ~isempty(aviObj)
            interpolatedTime = (0:1/videoFrameRate:max(timeVector))';
            
            % Interpolate the processed signal
            interpolatedSignal = interp1(timeVector, processedSignal, interpolatedTime, 'linear', 'extrap');
            
            % Check for finite values
            if any(~isfinite(interpolatedSignal))
                uialert(fig, 'Interpolated signal contains non-finite values.', 'Error');
                return;
            end
            
            % Update the signal plot with interpolated data
            cla(axSignal); % Clear existing plot
            if subtractLowPassEnabled
                plot(axSignal, interpolatedTime, interpolatedSignal, 'g', 'LineWidth', 1.5, 'DisplayName', 'High-Frequency Signal (Interpolated)');
            else
                plot(axSignal, interpolatedTime, interpolatedSignal, 'r', 'LineWidth', 1.5, 'DisplayName', 'Processed Signal (Interpolated)');
            end
            xlim(axSignal, [0, max(interpolatedTime)]);
            currentTimeLine = plot(axSignal, [0 0], ylim(axSignal), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Current Time');
            
            % Update the slider range based on the interpolated time
            timeSlider.Limits = [0, max(interpolatedTime)];
            timeSlider.Value = 0;  % Reset slider to start
            
            % Update the time label
            timeLabel.Text = sprintf('Time: %.2f s', timeSlider.Value);
            
            % Display the first frame of the video
            aviObj.CurrentTime = 0;
            try
                frame = readFrame(aviObj);
                imshow(frame, 'Parent', axVideo);
            catch ME
                uialert(fig, ['Failed to read the first frame of the video: ' ME.message], 'Error');
                return;
            end
        end
    end
    %% *** End of New Function ***
    
    %% New Callback Functions for Behavior Tagging
    
    function tagBehavior(behaviorType)
        % Get the current time from the slider
        currentTime = timeSlider.Value;
        
        % Append the current time to the appropriate behavior list
        switch behaviorType
            case 1
                behavior1Times = [behavior1Times; currentTime];
                uialert(fig, sprintf('Behavior 1 tagged at %.2f seconds.', currentTime), 'Behavior Tagged');
                
                % Plot a vertical line at the tagged time
                plot(axSignal, [currentTime, currentTime], ylim(axSignal), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Behavior 1');
                
            case 2
                behavior2Times = [behavior2Times; currentTime];
                uialert(fig, sprintf('Behavior 2 tagged at %.2f seconds.', currentTime), 'Behavior Tagged');
                
                % Plot a vertical line at the tagged time
                plot(axSignal, [currentTime, currentTime], ylim(axSignal), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Behavior 2');
        end
        
        % Update the legend to include new tags
        legend(axSignal, 'show');
        
        % Enable the "Save Tags" button if at least one tag exists
        if ~isempty(behavior1Times) || ~isempty(behavior2Times)
            saveTagsButton.Enable = 'on';
        end
    end

    function saveTags()
        % Ensure that a CSV file has been loaded
        if isempty(csvBaseName)
            uialert(fig, 'Please load a CSV file before saving tags.', 'Error');
            return;
        end
        
        % Ensure that there is at least one tag to save
        if isempty(behavior1Times) && isempty(behavior2Times)
            uialert(fig, 'No behavior tags to save.', 'Error');
            return;
        end
        
        % Define the Excel file name based on the CSV file name
        excelFileName = fullfile(csvPath, [csvBaseName, '_BehaviorTags.xlsx']);
        
        % Prepare data tables for each behavior
        if ~isempty(behavior1Times)
            behavior1Table = table(behavior1Times, 'VariableNames', {'Behavior1_Onset_Time_s'});
        else
            behavior1Table = table([], 'VariableNames', {'Behavior1_Onset_Time_s'});
        end
        
        if ~isempty(behavior2Times)
            behavior2Table = table(behavior2Times, 'VariableNames', {'Behavior2_Onset_Time_s'});
        else
            behavior2Table = table([], 'VariableNames', {'Behavior2_Onset_Time_s'});
        end
        
        % Write to Excel with separate sheets
        try
            if ~isempty(behavior1Times)
                writetable(behavior1Table, excelFileName, 'Sheet', 'Behavior1');
            end
            if ~isempty(behavior2Times)
                writetable(behavior2Table, excelFileName, 'Sheet', 'Behavior2');
            end
            uialert(fig, ['Behavior tags saved successfully to ', excelFileName], 'Success');
        catch ME
            uialert(fig, ['Failed to save tags: ' ME.message], 'Error');
            return;
        end
        
        % Now, save the axSignal figure as TIFF with overlaid tags
        try
            % Define the TIFF filename based on the CSV file name
            tiffFileName = fullfile(csvPath, [csvBaseName, '_Signal_With_Tags.tif']);
            
            % Create a temporary figure to hold the axSignal for saving
            tempFig = figure('Visible', 'off'); % Invisible figure
            copyobj(axSignal, tempFig); % Copy the axSignal to the temporary figure
            set(tempFig, 'Position', get(fig, 'Position')); % Match figure size
            
            % Adjust axes properties for better saving quality
            set(tempFig.Children, 'FontSize', 12); % Set font size
            title(tempFig.Children, get(axSignal.Title, 'String')); % Retain title
            xlabel(tempFig.Children, get(axSignal.XLabel, 'String'));
            ylabel(tempFig.Children, get(axSignal.YLabel, 'String'));
            
            % Save the temporary figure as TIFF
            exportgraphics(tempFig.Children, tiffFileName, 'Resolution', 300); % High resolution
            
            % Close the temporary figure
            close(tempFig);
            
            % Notify the user of successful save
            uialert(fig, ['Figure saved successfully as ', tiffFileName], 'Figure Saved');
        catch ME
            uialert(fig, ['Failed to save figure: ' ME.message], 'Error');
        end
    end
end

%% *** New Function: Apply Baseline Correction ***
function correctedSignal = applyBaselineCorrection(signal, timeVector)
    % This function removes a slow drift (baseline) from the signal
    % using cubic detrending. You can modify this function to use
    % more sophisticated baseline correction methods if needed.
    
    % Fit a cubic polynomial trend to the signal
    p = polyfit(timeVector, signal, 3); % Cubic fit for better trend removal
    
    % Calculate the baseline
    baseline = polyval(p, timeVector);
    
    % Subtract the baseline from the original signal
    correctedSignal = signal - baseline;
end
%% *** End of New Function ***
