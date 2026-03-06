function fiber_photometry_video_sync_v4_6
% Fiber Photometry + Video Sync (Single ROI)
% - Modern layout with uigridlayout
% - Zoomed-in view centered at current time (auto y-limit within window)
% - Baseline correction, LP smoothing, optional high-frequency (raw - LP)
% - Behavior tagging + export to Excel and TIFF snapshot
% - Peak detection with threshold = mean + 3*SD
% - Safe timer cleanup on close

clc; close all; warning off;

%% ===== Main figure & safe shutdown =====
fig = uifigure('Name','FP + Video Sync (Single ROI)', ...
    'Position',[60 60 1800 950]);
fig.CloseRequestFcn = @(~,~) cleanupAndClose();

%% ===== Top layout =====
root = uigridlayout(fig,[2 1]);
root.RowHeight   = {'1x', 120};
root.ColumnWidth = {'1x'};

% Content area
content = uigridlayout(root,[2 1]);
content.Layout.Row = 1;
content.Layout.Column = 1;
content.RowHeight   = {'2x','1x'};
content.ColumnWidth = {'1x'};
content.RowSpacing  = 12;
content.Padding     = [12 10 12 4];

% Video panel
videoCard = uipanel(content,'Title','Video','FontSize',13,'BorderType','line');
videoCard.Layout.Row = 1;
videoCard.Layout.Column = 1;

videoGrid = uigridlayout(videoCard,[1 1]);
videoGrid.Padding = [8 8 8 4];

axVideo = uiaxes(videoGrid);
axis(axVideo,'off');
axis(axVideo,'image');
axVideo.YDir = 'reverse';
axVideo.XDir = 'normal';

% Signal panel
signalCard = uipanel(content,'Title','470/410 Signal','FontSize',13,'BorderType','line');
signalCard.Layout.Row = 2;
signalCard.Layout.Column = 1;

signalGrid = uigridlayout(signalCard,[1 1]);
signalGrid.Padding = [8 8 8 4];

axSignal = uiaxes(signalGrid);
xlabel(axSignal,'Time (s)');
ylabel(axSignal,'Signal (a.u.)');
title(axSignal,'Processed Signal');
grid(axSignal,'on');
hold(axSignal,'on');

%% ===== Bottom control bar =====
ctrl = uipanel(root);
ctrl.Layout.Row = 2;
ctrl.Layout.Column = 1;

g = uigridlayout(ctrl,[2,22]);
g.RowHeight = {40, 44};
g.ColumnWidth   = { ...
    72,72,72,72,72,72, ...      % 1..6
    68,68,78,78, ...            % 7..10
    78,72,56,70, ...            % 11..14
    70,95,60,60, ...            % 15..18
    60,60,180,120};             % 19..22
g.Padding       = [8 6 8 6];
g.RowSpacing    = 4;
g.ColumnSpacing = 5;

% ----- Row 1 -----
timeSlider = uislider(g, ...
    'ValueChangingFcn', @(~,e) updateTime(e.Value), ...
    'ValueChangedFcn',  @(s,~) updateTime(s.Value));
timeSlider.Layout.Row = 1;
timeSlider.Layout.Column = [1 15];
timeSlider.Limits = [0 1];
timeSlider.Value  = 0;

peakBaseCB = uicheckbox(g,'Text','Custom baseline','Value',false, ...
    'Tooltip','Use a user-defined baseline time window for threshold calculation', ...
    'ValueChangedFcn', @(~,~) onPeakParamsChanged());
peakBaseCB.Layout.Row = 1;
peakBaseCB.Layout.Column = 16;

baseStartEdit = uieditfield(g,'numeric', ...
    'Limits',[0 Inf], ...
    'Value',0, ...
    'Enable','off', ...
    'Tooltip','Baseline start time (s)', ...
    'ValueChangedFcn', @(~,~) onPeakParamsChanged());
baseStartEdit.Layout.Row = 1;
baseStartEdit.Layout.Column = 17;

baseEndEdit = uieditfield(g,'numeric', ...
    'Limits',[0 Inf], ...
    'Value',10, ...
    'Enable','off', ...
    'Tooltip','Baseline end time (s)', ...
    'ValueChangedFcn', @(~,~) onPeakParamsChanged());
baseEndEdit.Layout.Row = 1;
baseEndEdit.Layout.Column = 18;

smoothLbl = uilabel(g,'Text','Smooth','HorizontalAlignment','right');
smoothLbl.Layout.Row = 1;
smoothLbl.Layout.Column = 19;

smoothEdit = uieditfield(g,'numeric', ...
    'Limits',[0 20], ...
    'Value',0.5, ...
    'Tooltip','Smoothing window (seconds) used only for peak detection', ...
    'ValueChangedFcn', @(~,~) onPeakParamsChanged());
smoothEdit.Layout.Row = 1;
smoothEdit.Layout.Column = 20;

timeLbl = uilabel(g,'Text','t = 0.00 s','HorizontalAlignment','right');
timeLbl.Layout.Row = 1;
timeLbl.Layout.Column = 21;

% ----- Row 2 -----
csvBtn = uibutton(g,'Text','Load CSV', ...
    'ButtonPushedFcn', @(~,~) loadCSV());
csvBtn.Layout.Row = 2;
csvBtn.Layout.Column = 1;

aviBtn = uibutton(g,'Text','Load Video', ...
    'Enable','off', ...
    'ButtonPushedFcn', @(~,~) loadAVI());
aviBtn.Layout.Row = 2;
aviBtn.Layout.Column = 2;

playBtn = uibutton(g,'Text','Play', ...
    'ButtonPushedFcn', @(~,~) playVideo());
playBtn.Layout.Row = 2;
playBtn.Layout.Column = 3;

pauseBtn = uibutton(g,'Text','Pause', ...
    'ButtonPushedFcn', @(~,~) pauseVideo());
pauseBtn.Layout.Row = 2;
pauseBtn.Layout.Column = 4;

backBtn = uibutton(g,'Text','<< 5s', ...
    'ButtonPushedFcn', @(~,~) moveTime(-5));
backBtn.Layout.Row = 2;
backBtn.Layout.Column = 5;

fwBtn = uibutton(g,'Text','5s >>', ...
    'ButtonPushedFcn', @(~,~) moveTime(5));
fwBtn.Layout.Row = 2;
fwBtn.Layout.Column = 6;

tagB1 = uibutton(g,'Text','Tag 1', ...
    'ButtonPushedFcn', @(~,~) tagBehavior(1));
tagB1.Layout.Row = 2;
tagB1.Layout.Column = 7;

tagB2 = uibutton(g,'Text','Tag 2', ...
    'ButtonPushedFcn', @(~,~) tagBehavior(2));
tagB2.Layout.Row = 2;
tagB2.Layout.Column = 8;

saveTagsButton = uibutton(g,'Text','Save', ...
    'ButtonPushedFcn', @(~,~) saveTags(), ...
    'Enable','on');
saveTagsButton.Layout.Row = 2;
saveTagsButton.Layout.Column = 9;

baselineCB = uicheckbox(g,'Text','Baseline','Value',false, ...
    'Tooltip','Cubic detrend on the low-pass signal', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
baselineCB.Layout.Row = 2;
baselineCB.Layout.Column = 10;

normCB = uicheckbox(g,'Text','Norm ΔF/F','Value',true, ...
    'Tooltip','Convert to ΔF/F% using robust F0 (10th percentile)', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
normCB.Layout.Row = 2;
normCB.Layout.Column = 11;

zoomCB = uicheckbox(g,'Text','Zoomed','Enable','off', ...
    'Tooltip','Show a window around the current time', ...
    'ValueChangedFcn', @(~,~) applyZoom(timeSlider.Value));
zoomCB.Layout.Row = 2;
zoomCB.Layout.Column = 12;

zoomSpin = uispinner(g, ...
    'Limits',[1 600], ...
    'Value',100, ...
    'Step',1, ...
    'Enable','off', ...
    'Tooltip','Zoom window width (s)', ...
    'ValueChangedFcn', @(~,~) applyZoom(timeSlider.Value));
zoomSpin.Layout.Row = 2;
zoomSpin.Layout.Column = 13;

lblCut = uilabel(g,'Text','LP Hz:','HorizontalAlignment','right');
lblCut.Layout.Row = 2;
lblCut.Layout.Column = 14;

cutEdit = uieditfield(g,'numeric', ...
    'Limits',[1e-3 100], ...
    'Value',5, ...
    'Tooltip','Butterworth low-pass cutoff (Hz)', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
cutEdit.Layout.Row = 2;
cutEdit.Layout.Column = 15;

subLPcb = uicheckbox(g,'Text','Subtract LP','Value',false, ...
    'Tooltip','High-frequency output = raw - low-pass', ...
    'ValueChangedFcn', @(~,~) onProcParamsChanged());
subLPcb.Layout.Row = 2;
subLPcb.Layout.Column = 16;

detectPeaksCB = uicheckbox(g,'Text','Detect peaks','Value',false, ...
    'Tooltip','Detect peaks using threshold = mean + 3*SD', ...
    'ValueChangedFcn', @(~,~) onPeakParamsChanged());
detectPeaksCB.Layout.Row = 2;
detectPeaksCB.Layout.Column = 17;

peakInfoLbl = uilabel(g,'Text','Peaks: off','HorizontalAlignment','left');
peakInfoLbl.Layout.Row = 2;
peakInfoLbl.Layout.Column = [18 22];

%% ===== Internal state =====
csvData = [];
timeVec = [];
sig470  = [];
sig410  = [];
sigRaw  = [];
sigLP   = [];
sigProc = [];

aviObj = [];
videoTimer = timer( ...
    'ExecutionMode','fixedRate', ...
    'Period',1/30, ...
    'TimerFcn',@(~,~) onTimer());

currentTimeLine = gobjects(1,1);
yFull = [];

csvBaseName = '';
csvPath = '';

behavior1Times = [];
behavior2Times = [];

peakTimes = [];
peakAmps = [];
peakAUCs = [];
peakFreqPerMin = NaN;
peakThreshold = NaN;
peakBaseMean = NaN;
peakBaseSD = NaN;
peakMarkers = gobjects(0);
peakThreshLine = gobjects(1,1);

%% ===== Callbacks =====

    function loadCSV()
        [f, p] = uigetfile('*.csv','Select the CSV file');
        if isequal(f,0), return; end

        csvPath = p;
        [~, csvBaseName, ~] = fileparts(f);

        try
            csvData = readtable(fullfile(p,f));
        catch ME
            uialert(fig,['Failed to read CSV: ' ME.message],'Error');
            return;
        end

        if width(csvData) < 3
            uialert(fig,'CSV must have at least 3 columns: Time, 410, 470.','Error');
            return;
        end

        TimeStamp = csvData{:,1};

        if iscell(TimeStamp) || isstring(TimeStamp) || ischar(TimeStamp)
            try
                TimeStamp = duration(string(TimeStamp), 'InputFormat','hh:mm:ss:SSSSSSS');
            catch
                try
                    TimeStamp = duration(string(TimeStamp));
                catch
                    uialert(fig,'Column 1 could not be parsed as duration/time.','Error');
                    return;
                end
            end
        end

        if ~isdatetime(TimeStamp) && ~isduration(TimeStamp)
            uialert(fig,'Column 1 must be datetime/duration or a parseable time string.','Error');
            return;
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
        if isempty(idx470), idx470 = min(3,width(csvData)); end

        sig410 = csvData{:,idx410};
        sig470 = csvData{:,idx470};

        sig410 = sig410(:);
        sig470 = sig470(:);

        finiteIdx = isfinite(timeVec) & isfinite(sig410) & isfinite(sig470) & sig410~=0;
        timeVec = timeVec(finiteIdx);
        sig410  = sig410(finiteIdx);
        sig470  = sig470(finiteIdx);

        sigRaw = sig470 ./ sig410;

        if isempty(timeVec)
            uialert(fig,'No valid data after cleaning.','Error');
            return;
        end

        processAndPlotSignal();

        timeSlider.Limits = [0 max(timeVec)];
        timeSlider.Value = 0;
        timeLbl.Text = sprintf('t = %.2f s',0);

        ensureTimeLine(0);

        aviBtn.Enable = 'on';
        zoomCB.Enable = 'on';
        zoomSpin.Enable = 'on';

        drawnow;
        yFull = ylim(axSignal);
    end

    function loadAVI()
        if isempty(timeVec) || isempty(sigProc)
            uialert(fig,'Load CSV first.','Error');
            return;
        end

        [f, p] = uigetfile({'*.avi;*.mp4;*.mov','Video Files'},'Select the video');
        if isequal(f,0), return; end

        try
            aviObj = VideoReader(fullfile(p,f));
        catch ME
            uialert(fig,['Failed to read video: ' ME.message],'Error');
            return;
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
        if isempty(aviObj)
            uialert(fig,'Load a video first.','Error');
            return;
        end
        if strcmp(videoTimer.Running,'off')
            start(videoTimer);
        end
    end

    function pauseVideo()
        if strcmp(videoTimer.Running,'on')
            stop(videoTimer);
        end
    end

    function onTimer()
        if isempty(aviObj), return; end

        if ~hasFrame(aviObj)
            stop(videoTimer);
            return;
        end

        try
            frame = readFrame(aviObj);
            imshow(frame,'Parent',axVideo);
        catch
            stop(videoTimer);
            return;
        end

        t = aviObj.CurrentTime;
        timeSlider.Value = min(max(t,timeSlider.Limits(1)), timeSlider.Limits(2));
        updateTime(timeSlider.Value);
    end

    function moveTime(dt)
        newT = min(max(timeSlider.Value + dt, timeSlider.Limits(1)), timeSlider.Limits(2));
        timeSlider.Value = newT;
        updateTime(newT);
    end

    function updateTime(t)
        if isempty(timeVec), return; end

        timeLbl.Text = sprintf('t = %.2f s',t);

        if ~isempty(aviObj)
            try
                aviObj.CurrentTime = min(max(t,0), max(0, aviObj.Duration - 1e-6));
                if hasFrame(aviObj)
                    frame = readFrame(aviObj);
                    imshow(frame,'Parent',axVideo);
                end
            catch
            end
        end

        ensureTimeLine(t);
        applyZoom(t);
    end

    function ensureTimeLine(t)
        yl = ylim(axSignal);
        if ~isgraphics(currentTimeLine)
            currentTimeLine = plot(axSignal,[t t],[yl(1) yl(2)],'k--','LineWidth',1.2);
        else
            currentTimeLine.XData = [t t];
            currentTimeLine.YData = [yl(1) yl(2)];
        end
    end

    function onProcParamsChanged()
        if isempty(timeVec), return; end
        processAndPlotSignal();
        ensureTimeLine(timeSlider.Value);
        applyZoom(timeSlider.Value);
    end

    function processAndPlotSignal()
        dt = median(diff(timeVec));
        Fs = 1 / max(dt, eps);

        fc = cutEdit.Value;
        if fc >= Fs/2
            fc = max(1e-3, Fs/2 - 1e-3);
            cutEdit.Value = fc;
        end

        [b,a] = butter(4, fc/(Fs/2), 'low');
        sigLP = filtfilt(b,a,sigRaw);

        if baselineCB.Value
            sigLP = applyBaselineCorrection(sigLP, timeVec);
        end

        if subLPcb.Value
            sigProc = sigRaw - sigLP;
        else
            sigProc = sigLP;
        end

        if normCB.Value
            base = prctile(sigProc,10);
            if ~isfinite(base) || abs(base) < eps
                base = eps;
            end
            sigProc = (sigProc - base) ./ base * 100;
        end

        cla(axSignal);
        plot(axSignal,timeVec,sigProc,'LineWidth',1.2);
        xlim(axSignal,[0 max(timeVec)]);

        if subLPcb.Value
            title(axSignal,'High-Frequency (Raw - LP)');
        else
            title(axSignal,'Processed Signal');
        end

        xlabel(axSignal,'Time (s)');
        if normCB.Value
            ylabel(axSignal,'\DeltaF/F (%)');
        else
            ylabel(axSignal,'470/410 (a.u.)');
        end
        grid(axSignal,'on');

        drawnow;
        yFull = ylim(axSignal);

        detectAndPlotPeaks();
        ensureTimeLine(timeSlider.Value);
    end

    function applyZoom(tCenter)
        if isempty(timeVec) || ~isvalid(zoomCB) || ~zoomCB.Value
            xlim(axSignal,[timeSlider.Limits(1) timeSlider.Limits(2)]);
            if ~isempty(sigProc)
                ymin = min(sigProc);
                ymax = max(sigProc);
                if ~isfinite(ymin) || ~isfinite(ymax) || ymin == ymax
                    pad = 0.5;
                    ymin = ymin - pad;
                    ymax = ymax + pad;
                end
                ylim(axSignal,[ymin ymax]);
                yFull = [ymin ymax];
            elseif ~isempty(yFull)
                ylim(axSignal,yFull);
            end
            ensureTimeLine(timeSlider.Value);
            return;
        end

        W = max(1, zoomSpin.Value);
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

        idx = (timeVec >= t1 & timeVec <= t2);
        if any(idx)
            ymin = min(sigProc(idx));
            ymax = max(sigProc(idx));
            if ~isfinite(ymin) || ~isfinite(ymax) || ymin == ymax
                pad = 0.5;
                ymin = ymin - pad;
                ymax = ymax + pad;
            end
            ylim(axSignal,[ymin ymax]);
            if isgraphics(currentTimeLine)
                currentTimeLine.YData = [ymin ymax];
            end
        end
    end

    function onPeakParamsChanged()
        if peakBaseCB.Value
            baseStartEdit.Enable = 'on';
            baseEndEdit.Enable   = 'on';
        else
            baseStartEdit.Enable = 'off';
            baseEndEdit.Enable   = 'off';
        end

        if isempty(timeVec) || isempty(sigProc)
            return;
        end

        detectAndPlotPeaks();
        ensureTimeLine(timeSlider.Value);
        applyZoom(timeSlider.Value);
    end

    function detectAndPlotPeaks()
        try
            if ~isempty(peakMarkers)
                delete(peakMarkers(isgraphics(peakMarkers)));
            end
        catch
        end
        peakMarkers = gobjects(0);

        try
            if isgraphics(peakThreshLine)
                delete(peakThreshLine);
            end
        catch
        end

        peakTimes = [];
        peakAmps = [];
        peakAUCs = [];
        peakFreqPerMin = NaN;
        peakThreshold = NaN;
        peakBaseMean = NaN;
        peakBaseSD = NaN;

        if isempty(timeVec) || isempty(sigProc)
            peakInfoLbl.Text = 'Peaks: no data';
            return;
        end

        if ~detectPeaksCB.Value
            peakInfoLbl.Text = 'Peaks: off';
            return;
        end

        y = sigProc(:);
        t = timeVec(:);

        dt = median(diff(t));
        if ~isfinite(dt) || dt <= 0
            dt = 0.05;
        end

        smoothSec = smoothEdit.Value;
        winN = max(1, round(smoothSec / dt));

        if winN > 1
            yDet = smoothdata(y,'movmean',winN);
        else
            yDet = y;
        end

        if peakBaseCB.Value
            t1 = baseStartEdit.Value;
            t2 = baseEndEdit.Value;

            if t2 <= t1
                peakInfoLbl.Text = 'Peaks: invalid baseline';
                return;
            end

            idxBase = (t >= t1 & t <= t2);
            if ~any(idxBase)
                peakInfoLbl.Text = 'Peaks: empty baseline';
                return;
            end
        else
            idxBase = true(size(t));
        end

        baseData = yDet(idxBase);
        peakBaseMean = mean(baseData,'omitnan');
        peakBaseSD   = std(baseData,0,'omitnan');
        peakThreshold = peakBaseMean + 2*peakBaseSD;

        if ~isfinite(peakThreshold)
            peakInfoLbl.Text = 'Peaks: bad threshold';
            return;
        end

        minPeakDistSec = max(0.2, smoothSec);

        [pks, locs] = findpeaks(yDet, t, ...
            'MinPeakHeight', peakThreshold, ...
            'MinPeakDistance', minPeakDistSec);

        hold(axSignal,'on');
        peakThreshLine = yline(axSignal, peakThreshold, 'm--', 'LineWidth', 1.1);

        if isempty(pks)
            peakInfoLbl.Text = sprintf('Peaks:0 | Thr:%.3f', peakThreshold);
            return;
        end

        peakIdx = zeros(size(locs));
        for i = 1:numel(locs)
            [~, peakIdx(i)] = min(abs(t - locs(i)));
        end

        nPeaks = numel(peakIdx);
        peakTimes = t(peakIdx);
        peakAmps  = nan(nPeaks,1);
        peakAUCs  = nan(nPeaks,1);

        for i = 1:nPeaks
            ip = peakIdx(i);

            iL = ip;
            while iL > 1 && yDet(iL) > peakThreshold
                iL = iL - 1;
            end

            iR = ip;
            while iR < numel(yDet) && yDet(iR) > peakThreshold
                iR = iR + 1;
            end

            peakAmps(i) = yDet(ip) - peakBaseMean;

            segT = t(iL:iR);
            segY = yDet(iL:iR) - peakBaseMean;
            segY(segY < 0) = 0;
            peakAUCs(i) = trapz(segT, segY);
        end

        durationMin = (t(end) - t(1)) / 60;
        if durationMin > 0
            peakFreqPerMin = nPeaks / durationMin;
        else
            peakFreqPerMin = NaN;
        end

        peakMarkers = plot(axSignal, peakTimes, y(peakIdx), 'mo', ...
            'MarkerSize', 5, 'LineWidth', 1.0);

        peakInfoLbl.Text = sprintf( ...
            'Peaks:%d | Freq:%.2f/min | Amp:%.3f | AUC:%.3f', ...
            nPeaks, peakFreqPerMin, ...
            mean(peakAmps,'omitnan'), ...
            mean(peakAUCs,'omitnan'));
    end

    function tagBehavior(whichOne)
        t = timeSlider.Value;
        switch whichOne
            case 1
                behavior1Times = [behavior1Times; t];
                plot(axSignal,[t t],ylim(axSignal),'r--','LineWidth',1.1);
                uialert(fig,sprintf('Behavior 1 tagged at %.2f s',t),'Tagged');
            case 2
                behavior2Times = [behavior2Times; t];
                plot(axSignal,[t t],ylim(axSignal),'b--','LineWidth',1.1);
                uialert(fig,sprintf('Behavior 2 tagged at %.2f s',t),'Tagged');
        end
        saveTagsButton.Enable = 'on';
    end

    function saveTags()
        if isempty(csvBaseName)
            uialert(fig,'Load a CSV before saving.','Error');
            return;
        end

        if isempty(behavior1Times) && isempty(behavior2Times) && isempty(peakTimes)
            uialert(fig,'No tags or peaks to save.','Error');
            return;
        end

        xlsx = fullfile(csvPath,[csvBaseName '_BehaviorTags.xlsx']);

        try
            if ~isempty(behavior1Times)
                writetable(table(behavior1Times,'VariableNames',{'Behavior1_Onset_Time_s'}), ...
                    xlsx,'Sheet','Behavior1');
            end

            if ~isempty(behavior2Times)
                writetable(table(behavior2Times,'VariableNames',{'Behavior2_Onset_Time_s'}), ...
                    xlsx,'Sheet','Behavior2');
            end

            if ~isempty(peakTimes)
                peakTbl = table(peakTimes, peakAmps, peakAUCs, ...
                    'VariableNames', {'PeakTime_s','PeakAmplitude','PeakAUC'});
                writetable(peakTbl, xlsx, 'Sheet', 'Peaks');

                summaryTbl = table(peakFreqPerMin, peakBaseMean, peakBaseSD, peakThreshold, ...
                    mean(peakAmps,'omitnan'), mean(peakAUCs,'omitnan'), ...
                    'VariableNames', {'PeakFreq_per_min','BaselineMean','BaselineSD', ...
                                      'Threshold','MeanPeakAmplitude','MeanPeakAUC'});
                writetable(summaryTbl, xlsx, 'Sheet', 'PeakSummary');
            end
        catch ME
            uialert(fig,['Failed to write Excel: ' ME.message],'Error');
            return;
        end

        try
            tiffName = fullfile(csvPath,[csvBaseName '_Signal_With_Tags.tif']);
            exportgraphics(axSignal, tiffName, 'Resolution', 300);
            uialert(fig, sprintf('Saved:\n%s\n%s', xlsx, tiffName), 'Success');
        catch ME
            uialert(fig,['Failed to save TIFF: ' ME.message],'Error');
        end
    end

%% ===== Helpers =====
    function corrected = applyBaselineCorrection(signal, tvec)
        p = polyfit(tvec, signal, 3);
        corrected = signal - polyval(p, tvec);
    end

    function cleanupAndClose()
        try
            if ~isempty(videoTimer) && isvalid(videoTimer)
                stop(videoTimer);
                delete(videoTimer);
            end
        catch
        end
        delete(fig);
    end

end