function fiber_photometry_video_sync_v4_4
% Fiber Photometry + Video Sync (Single ROI)
% - Modern layout with uigridlayout
% - Zoomed-in view centered at current time (auto y-limit within window)
% - Baseline correction, LP smoothing, optional high-frequency (raw - LP)
% - Behavior tagging + export to Excel (two sheets) and TIFF snapshot
% - Safe timer cleanup on close

clc; close all; warning off;

%% ===== Main figure & safe shutdown =====
fig = uifigure('Name','FP + Video Sync (Single ROI)','Position',[60 60 1500 950]);
fig.CloseRequestFcn = @(~,~) cleanupAndClose();

%% ===== Top layout: left (signal), right (video) =====
% Use two main rows:
%   Row 1 = content (video + traces)
%   Row 2 = controls (slider + buttons)
root = uigridlayout(fig,[2 1]);
root.RowHeight   = {'1x', 120};   % Main content, then controls
root.ColumnWidth = {'1x'};

% Content area: 2 rows (video 2/3, signal 1/3)
content = uigridlayout(root,[2 1]);
content.Layout.Row = 1; content.Layout.Column = 1;
content.RowHeight   = {'2x','1x'};   % Top 2/3 for video, bottom 1/3 for traces
content.ColumnWidth = {'1x'};
content.RowSpacing  = 12;
content.Padding     = [12 10 12 4];

% Top: Video card
videoCard = uipanel(content,'Title','Video','FontSize',13,'BorderType','line');
videoCard.Layout.Row = 1; videoCard.Layout.Column = 1;
videoGrid = uigridlayout(videoCard,[1 1]); videoGrid.Padding=[8 8 8 4];
axVideo = uiaxes(videoGrid); axis(axVideo,'off'); axis(axVideo,'image');
axVideo.YDir = 'reverse'; axVideo.XDir = 'normal';

% Bottom: Signal card
signalCard = uipanel(content,'Title','470/410 Signal','FontSize',13,'BorderType','line');
signalCard.Layout.Row = 2; signalCard.Layout.Column = 1;
signalGrid = uigridlayout(signalCard,[1 1]); signalGrid.Padding=[8 8 8 4];
axSignal = uiaxes(signalGrid);
xlabel(axSignal,'Time (s)'); ylabel(axSignal,'Signal (a.u.)');
title(axSignal,'Processed Signal'); grid(axSignal,'on'); hold(axSignal,'on');


%% ===== Bottom control bar (2 rows) =====
ctrl = uipanel(root); 
ctrl.Layout.Row = 2; 
ctrl.Layout.Column = 1;

% 20 columns: Row 1 = thin slider row, Row 2 = controls
g = uigridlayout(ctrl,[2,16]);
g.RowHeight     = {40, 50};      % Row 1 (slider+compact controls), Row 2 (buttons)
g.ColumnWidth   = { ...
    80,80,80,80,80,80, ...       % 1..6   (left cluster / slider)
    80, ...                    % 7      (flex for slider stretch)
    80,80,80, ...                % 8..10  (still under slider)
    80, ...                      % 11     (slider end)
    80,80,80, ...             % 12..14 Normalize, Zoomed in, Zoom width
    80,80, ...                 % 15..16 LP label, LP value
    };                        % 20     Time label (a bit wider looks nicer)
g.Padding       = [10 8 10 8];
g.RowSpacing    = 6;
g.ColumnSpacing = 8;

% ---------- Row 1: Slider (1..11) + compact controls (12..19) + time label (20)
timeSlider  = uislider(g,'ValueChangingFcn',@(s,e) updateTime(e.Value), ...
                          'ValueChangedFcn', @(s,~) updateTime(s.Value));
timeSlider.Layout.Row = 1;
timeSlider.Layout.Column = [1 15];
% 20) Time label at end of Row 1
timeLbl = uilabel(g,'Text','Time: 0.00 s','HorizontalAlignment','right');
timeLbl.Layout.Row = 1; 
timeLbl.Layout.Column = 16;





% ---------- Row 2: Files/Playback (left cluster)
csvBtn   = uibutton(g,'Text','Load CSV','ButtonPushedFcn', @(~,~) loadCSV()); 
csvBtn.Layout.Row = 2; 
csvBtn.Layout.Column = 1;

aviBtn   = uibutton(g,'Text','Load Video','Enable','off','ButtonPushedFcn', @(~,~) loadAVI());
aviBtn.Layout.Row = 2; 
aviBtn.Layout.Column = 2;

playBtn  = uibutton(g,'Text','Play','ButtonPushedFcn', @(~,~) playVideo());
playBtn.Layout.Row = 2; 
playBtn.Layout.Column = 3;

pauseBtn = uibutton(g,'Text','Pause','ButtonPushedFcn', @(~,~) pauseVideo());
pauseBtn.Layout.Row = 2; 
pauseBtn.Layout.Column = 4;

backBtn  = uibutton(g,'Text','<< 5s','ButtonPushedFcn', @(~,~) moveTime(-5));
backBtn.Layout.Row = 2; 
backBtn.Layout.Column = 5;

fwBtn    = uibutton(g,'Text','5s >>','ButtonPushedFcn', @(~,~) moveTime(5));
fwBtn.Layout.Row  = 2;  
fwBtn.Layout.Column = 6;

% Tag buttons (Row 2)
tagB1 = uibutton(g,'Text','Tag Behavior 1','ButtonPushedFcn', @(~,~) tagBehavior(1));
tagB1.Layout.Row = 2; 
tagB1.Layout.Column = 7;

tagB2 = uibutton(g,'Text','Tag Behavior 2','ButtonPushedFcn', @(~,~) tagBehavior(2));
tagB2.Layout.Row = 2; 
tagB2.Layout.Column = 8;

% Save Tags (use the grid, not absolute position)
saveTagsButton = uibutton(g,'Text','Save Tags', ...
    'ButtonPushedFcn', @(~,~) saveTags(), 'Enable','on');
saveTagsButton.Layout.Row = 2; 
saveTagsButton.Layout.Column = 9;

% Baseline correction toggle
baselineCB = uicheckbox(g,'Text','Baseline Corr','Value',false, ...
    'Tooltip','Cubic detrend on the LP signal', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
baselineCB.Layout.Row = 2; 
baselineCB.Layout.Column = 10;


% 12) Normalize (ΔF/F%)
normCB = uicheckbox(g,'Text','Normalize (ΔF/F%)','Value',true, ...
    'Tooltip','Convert to ΔF/F% using robust F0 (10th percentile)', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
normCB.Layout.Row = 2; 
normCB.Layout.Column = 11;

% 13) Zoomed in
zoomCB = uicheckbox(g,'Text','Zoomed in','Enable','off', ...
    'Tooltip','Show a window around the current time', ...
    'ValueChangedFcn', @(~,~) applyZoom(timeSlider.Value));
zoomCB.Layout.Row = 2; 
zoomCB.Layout.Column = 12;

% 14) Zoom Width (s)
zoomSpin = uispinner(g,'Limits',[1 600],'Value',100,'Step',1,'Enable','off', ...
    'Tooltip','Zoom window width (seconds)', ...
    'ValueChangedFcn', @(~,~) applyZoom(timeSlider.Value));
zoomSpin.Layout.Row = 2; 
zoomSpin.Layout.Column = 13;

% 15) LP Cutoff label
lblCut = uilabel(g,'Text','LP Cutoff (Hz):','HorizontalAlignment','right');
lblCut.Layout.Row = 2; 
lblCut.Layout.Column = 14;

% 16) LP Cutoff value
cutEdit = uieditfield(g,'numeric','Limits',[1e-3 100],'Value',5, ...
    'Tooltip','Butterworth low-pass cutoff (Hz)', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
cutEdit.Layout.Row = 2; 
cutEdit.Layout.Column = 15;

% 17) Subtract Low-Pass
subLPcb = uicheckbox(g,'Text','Subtract Low-Pass','Value',false, ...
    'Tooltip','Output high-frequency = raw - low-pass', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
subLPcb.Layout.Row = 2; 
subLPcb.Layout.Column = 16;



%% ===== Internal state =====
csvData = [];                   % raw table
timeVec = [];                   % seconds, column vector
sig470  = []; sig410 = [];      % raw channels (optional to keep)
sigRaw  = [];                   % 470/410
sigLP   = [];                   % low-pass
sigProc = [];                   % final processed (LP or Raw-LP)
aviObj  = [];                   % VideoReader
videoTimer = timer('ExecutionMode','fixedRate','Period',1/30,'TimerFcn',@(~,~) onTimer());
currentTimeLine = gobjects(1,1);% vertical line handle
yFull = [];                     % cached y-lims for full data
csvBaseName = ''; csvPath = ''; % for saving outputs
behavior1Times = []; behavior2Times = [];

%% ===== Callbacks =====

    function loadCSV()
        [f, p] = uigetfile('*.csv','Select the CSV file');
        if isequal(f,0), return; end
        csvPath = p;
        [~, csvBaseName, ~] = fileparts(f);

        try
            csvData = readtable(fullfile(p,f));
        catch ME
            uialert(fig,['Failed to read CSV: ' ME.message],'Error'); return;
        end

        % Expect time as datetime/duration in col 1; 2=410, 3/4=470
        if width(csvData) < 3
            uialert(fig,'CSV must have ≥3 columns (Time, 410, 470).','Error'); return;
        end

        TimeStamp = csvData{:,1};
        if ~isdatetime(TimeStamp) && ~isduration(TimeStamp)
            uialert(fig,'Column 1 must be datetime/duration.','Error'); return;
        end
        if isduration(TimeStamp)
            timeVec = seconds(TimeStamp) - seconds(TimeStamp(1));
        else
            timeVec = seconds(TimeStamp - TimeStamp(1));
        end
        timeVec = timeVec(:);

        colNames = lower(string(csvData.Properties.VariableNames));
        idx410 = find(contains(colNames,"410"),1,'first');
        idx470 = find(contains(colNames,"470"),1,'first');
        if isempty(idx410), idx410 = 2; end
        if isempty(idx470), idx470 = min(4,width(csvData)); end

        sig410 = csvData{:,idx410};
        sig470 = csvData{:,idx470};
        sigRaw = (sig470(:) ./ sig410(:));

        finiteIdx = isfinite(timeVec) & isfinite(sigRaw);
        timeVec = timeVec(finiteIdx);
        sigRaw  = sigRaw(finiteIdx);

        if isempty(timeVec)
            uialert(fig,'No valid data after cleaning.','Error'); return;
        end

        % Process & plot
        processAndPlotSignal();

        % Slider range
        timeSlider.Limits = [0 max(timeVec)];
        timeSlider.Value  = 0;
        timeLbl.Text = sprintf('Time: %.2f s',0);

        % current time line
        ensureTimeLine(0);

        % Enable video button; enable zoom toggles (only after data loaded)
        aviBtn.Enable = 'on';
        zoomCB.Enable = 'on';
        zoomSpin.Enable = 'on';

        % Cache full-range y for quick restore
        drawnow;
        yFull = ylim(axSignal);
    end

    function loadAVI()
        if isempty(timeVec) || isempty(sigProc)
            uialert(fig,'Load CSV first.','Error'); return;
        end
        [f, p] = uigetfile({'*.avi;*.mp4;*.mov','Video Files'},'Select the video');
        if isequal(f,0), return; end
        try
            aviObj = VideoReader(fullfile(p,f));
        catch ME
            uialert(fig,['Failed to read video: ' ME.message],'Error'); return;
        end
        videoTimer.Period = 1 / aviObj.FrameRate;
        aviObj.CurrentTime = 0;
        try
            frame = readFrame(aviObj);
            imshow(frame,'Parent',axVideo);
        catch ME
            uialert(fig,['Cannot read first frame: ' ME.message],'Error');
        end
    end

    function playVideo()
        if isempty(aviObj), uialert(fig,'Load a video first.','Error'); return; end
        if strcmp(videoTimer.Running,'off'), start(videoTimer); end
    end

    function pauseVideo()
        if strcmp(videoTimer.Running,'on'), stop(videoTimer); end
    end

    function onTimer()
        if isempty(aviObj), return; end
        if ~hasFrame(aviObj)
            stop(videoTimer); return;
        end
        try
            frame = readFrame(aviObj);
            imshow(frame,'Parent',axVideo);
        catch
            stop(videoTimer); return;
        end
        t = aviObj.CurrentTime;
        % update slider/plot (throttled by timer period)
        timeSlider.Value = t;
        updateTime(t);
    end

    function moveTime(dt)
        newT = min(max(timeSlider.Value + dt, timeSlider.Limits(1)), timeSlider.Limits(2));
        timeSlider.Value = newT;
        updateTime(newT);
    end

    function updateTime(t)
        if isempty(timeVec), return; end
        timeLbl.Text = sprintf('Time: %.2f s',t);
        % update video frame if video loaded
        if ~isempty(aviObj)
            try
                aviObj.CurrentTime = min(max(t,0), aviObj.Duration - 1e-6);
                if hasFrame(aviObj)
                    frame = readFrame(aviObj);
                    imshow(frame,'Parent',axVideo);
                end
            catch
                % ignore
            end
        end
        % update vertical line & zoom window
        ensureTimeLine(t);
        applyZoom(t);
    end

    function ensureTimeLine(t)
        if ~isgraphics(currentTimeLine)
            currentTimeLine = plot(axSignal,[t t], ylim(axSignal), 'k--','LineWidth',1.2);
        else
            currentTimeLine.XData = [t t];
            % Refresh Y to current axes limits
            yl = ylim(axSignal);
            currentTimeLine.YData = [yl(1) yl(2)];
        end
    end

    function onProcParamsChanged()
        if isempty(timeVec), return; end
        processAndPlotSignal();
        % keep current time line & zoom consistent
        ensureTimeLine(timeSlider.Value);
        applyZoom(timeSlider.Value);
    end

    function processAndPlotSignal()
        % Determine Fs and validate LP cutoff
        dt = median(diff(timeVec));
        Fs = 1/max(dt, eps);
        fc = cutEdit.Value;
        if fc >= Fs/2
            fc = max(1e-3, Fs/2 - 1e-3);
            cutEdit.Value = fc;
        end

        % Low-pass (zero-phase)
        [b,a] = butter(4, fc/(Fs/2), 'low');
        sigLP = filtfilt(b,a,sigRaw);

        % Baseline correction (on the LP result)
        if baselineCB.Value
            sigLP = applyBaselineCorrection(sigLP, timeVec);
        end

        % Choose final signal
        if subLPcb.Value
            sigProc = sigRaw - sigLP;     % high-frequency component
        else
            sigProc = sigLP;              % smoothed (and baseline corrected)
        end

        if normCB.Value
            base = prctile(sigProc, 10);              % robust baseline
            if base == 0 || ~isfinite(base), base = eps; end
            sigProc = (sigProc - base) ./ base * 100; % ΔF/F%
        end

        % Draw
        cla(axSignal);
        plot(axSignal, timeVec, sigProc, 'LineWidth',1.4);
        xlim(axSignal,[0 max(timeVec)]);
        % Title
        if subLPcb.Value
            ttl = 'High-Frequency (Raw - LP)';
        else
            ttl = 'Processed Signal';
        end
        title(axSignal, ttl);

        % Axis labels
        xlabel(axSignal,'Time (s)');
        if normCB.Value
            ylabel(axSignal,'\DeltaF/F (%)');
        else
            ylabel(axSignal,'470/410 (a.u.)');
        end

        grid(axSignal,'on');

        % Recreate current-time line on top
        ensureTimeLine(timeSlider.Value);

        % Update cached full y-lims
        drawnow;
        yFull = ylim(axSignal);
    end

    function applyZoom(tCenter)
        if isempty(timeVec) || ~isvalid(zoomCB) || ~zoomCB.Value
            % ---- Zoom OFF: full X and global Y from current processed signal ----
            xlim(axSignal, [timeSlider.Limits(1) timeSlider.Limits(2)]);
            if ~isempty(sigProc)
                % recompute robust global Y from current signal (safer than stale yFull)
                ymin = min(sigProc); ymax = max(sigProc);
                if ~isfinite(ymin) || ~isfinite(ymax) || ymin==ymax
                    pad = 0.5; ymin = ymin - pad; ymax = ymax + pad;
                end
                ylim(axSignal, [ymin ymax]);
                yFull = [ymin ymax];  % refresh the cache
            elseif ~isempty(yFull)
                ylim(axSignal, yFull);
            end
            ensureTimeLine(timeSlider.Value);
            return;
        end

        % ---- Zoom ON: set [t1,t2] and auto-Y to local min/max in window ----
        W  = max(1, zoomSpin.Value);
        hw = W/2;
        t1 = max(timeSlider.Limits(1), tCenter - hw);
        t2 = min(timeSlider.Limits(2), tCenter + hw);
        if (t2 - t1) < W
            if t1 <= timeSlider.Limits(1)
                t2 = min(timeSlider.Limits(2), t1 + W);
            elseif t2 >= timeSlider.Limits(2)
                t1 = max(timeSlider.Limits(1), t2 - W);
            end
        end
        xlim(axSignal,[t1 t2]);

        % Auto y within window
        idx = (timeVec >= t1 & timeVec <= t2);
        if any(idx)
            ymin = min(sigProc(idx)); ymax = max(sigProc(idx));
            if ~isfinite(ymin) || ~isfinite(ymax) || ymin==ymax
                pad = 0.5; ymin = ymin - pad; ymax = ymax + pad;
            end
            ylim(axSignal,[ymin ymax]);
            if isgraphics(currentTimeLine)
                currentTimeLine.YData = [ymin ymax];
            end
        end
    end

    function tagBehavior(whichOne)
        t = timeSlider.Value;
        switch whichOne
            case 1
                behavior1Times = [behavior1Times; t];
                plot(axSignal,[t t], ylim(axSignal), 'r--','LineWidth',1.2);
                uialert(fig,sprintf('Behavior 1 tagged at %.2f s',t),'Tagged');
            case 2
                behavior2Times = [behavior2Times; t];
                plot(axSignal,[t t], ylim(axSignal), 'b--','LineWidth',1.2);
                uialert(fig,sprintf('Behavior 2 tagged at %.2f s',t),'Tagged');
        end
        saveBtn.Enable = 'on';
    end

    function saveTags()
        if isempty(csvBaseName)
            uialert(fig,'Load a CSV before saving.','Error'); return;
        end
        if isempty(behavior1Times) && isempty(behavior2Times)
            uialert(fig,'No tags to save.','Error'); return;
        end

        xlsx = fullfile(csvPath, [csvBaseName '_BehaviorTags.xlsx']);
        try
            if ~isempty(behavior1Times)
                writetable(table(behavior1Times,'VariableNames',{'Behavior1_Onset_Time_s'}), ...
                    xlsx,'Sheet','Behavior1');
            end
            if ~isempty(behavior2Times)
                writetable(table(behavior2Times,'VariableNames',{'Behavior2_Onset_Time_s'}), ...
                    xlsx,'Sheet','Behavior2');
            end
        catch ME
            uialert(fig,['Failed to write Excel: ' ME.message],'Error'); return;
        end

        % Save a TIFF of current axSignal
        try
            tiffName = fullfile(csvPath,[csvBaseName '_Signal_With_Tags.tif']);
            tmp = figure('Visible','off'); axTmp = copyobj(axSignal,tmp);
            set(tmp,'Position',[100 100 1200 420]);
            exportgraphics(axTmp, tiffName, 'Resolution', 300);
            close(tmp);
            uialert(fig, sprintf('Saved:\n%s\n%s', xlsx, tiffName), 'Success');
        catch ME
            uialert(fig,['Failed to save TIFF: ' ME.message],'Error');
        end
    end

%% ===== Helpers =====
    function corrected = applyBaselineCorrection(signal, tvec)
        % cubic detrend
        p = polyfit(tvec, signal, 3);
        corrected = signal - polyval(p, tvec);
    end

    function cleanupAndClose()
        try
            if isvalid(videoTimer)
                stop(videoTimer); delete(videoTimer);
            end
        catch
        end
        delete(fig);
    end

end

