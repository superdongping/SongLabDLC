function fiber_photometry_video_sync_v4
    % Create a UI figure
    fig = uifigure('Name', 'Fiber Photometry and Video Synchronization', 'Position', [50, 50, 1200, 1000]); % Increased width to accommodate new buttons
    
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
    
    % Add a slider to scroll through time
    timeSlider = uislider(fig, 'Position', [60, 150, 915, 20], ... % Adjusted width to match increased figure size
        'ValueChangedFcn', @(sld, event) updateTime());
    timeSlider.Limits = [0, 1];  % Initialize with 0-1, will be updated after loading data
    
    % Display current time label
    timeLabel = uilabel(fig, 'Position', [400, 450, 200, 20], 'Text', 'Time: 0.00 s', ...
        'FontSize', 12, 'HorizontalAlignment', 'center'); % Adjusted position for larger figure
    
    % Axes for displaying the video
    axVideo = uiaxes(fig, 'Position', [100, 420, 900, 600]);
    axis(axVideo, 'off');  % Hide axes for better video display
    
    % Axes for plotting 470/410 signal
    axSignal = uiaxes(fig, 'Position', [20, 180, 960, 250]);
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
        requiredColumns = [1, 2, 4];
        if width(csvData) < max(requiredColumns)
            uialert(fig, 'CSV file does not contain the required columns (at least 4 columns expected).', 'Error');
            return;
        end
        
        % Extract timestamps and LED signals
        try
            TimeStamp = csvData{:, 1};  % Extracted timestamps as datetime or duration type
            LED_410 = csvData{:, 2};
            LED_470 = csvData{:, 4};
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

        % Apply Low-Pass Butterworth Filter to Smooth the Signal
        Fc = 5;             % Cutoff frequency in Hz (adjust as needed), change to 0.1 smooth the signal
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

        % Plot the smoothed signal on axSignal
        cla(axSignal); % Clear existing plot
        plot(axSignal, timeVector, Real_Signal_Smoothed, 'r', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Signal');
        
        % Set x-axis limits to match slider limits
        xlim(axSignal, [0, max(timeVector)]);
        
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
        try
            interpolatedSignal = interp1(timeVector, Real_Signal_Smoothed, interpolatedTime, 'linear', 'extrap');
        catch ME
            uialert(fig, ['Interpolation failed: ' ME.message], 'Error');
            return;
        end

        % Check for finite values in interpolatedSignal
        if any(~isfinite(interpolatedSignal))
            uialert(fig, 'Interpolated signal contains non-finite values.', 'Error');
            return;
        end

        % Update the signal plot with interpolated data
        cla(axSignal); % Clear existing plot
        plot(axSignal, interpolatedTime, interpolatedSignal, 'r', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Signal');
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
