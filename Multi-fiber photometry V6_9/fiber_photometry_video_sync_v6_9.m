function fiber_photometry_video_sync_v6_9

% GUI for Fiber Photometry Multi-ROI with Real-Time Corr & Lag (optimized, robust)
% - Clean bottom layout (slider row + controls row)
% - Persistent plot handles; update X/YData only (no replot while playing)
% - Cached preprocessing (Smooth, ΔF/F with robust baseline)
% - Upper-triangle xcorr; throttled by wall-clock
% - Safe cleanup & defensive handle checks to avoid "Invalid or deleted object"

clc; warning off; close all;

%% Main UI figure
fig = uifigure('Name','FP + Video Sync (Real-Time Corr/Lag)','Position',[50,50,1500,950]);
fig.CloseRequestFcn = @(~,~) cleanupAndClose();  % safe shutdown

%% Panels for traces (left), video (right-top), corr/lag (right-bottom)
panelTr = uipanel(fig,'Position',[20,140,740,780]);
gridTr  = uigridlayout(panelTr,[8,1],'RowHeight',repmat({'1x'},1,8));
axTr    = gobjects(1,8);
for i = 1:8
    axTr(i) = uiaxes(gridTr);
    axTr(i).XTick = []; axTr(i).XColor = 'none';
    title(axTr(i),sprintf('ROI %d',i),'FontSize',12);
end
xlabel(axTr(8),'Time (s)'); axTr(8).XColor = 'k'; axTr(8).XTickMode = 'auto';

axVideo = uiaxes(fig,'Position',[780,520,700,400]); axis(axVideo,'off');
axVideo.YDir = 'reverse';    % top-left of the frame appears at top-left in the axes
axVideo.XDir = 'normal';     % no left-right mirroring
axis(axVideo,'image');       % correct aspect ratio

panelCL = uipanel(fig,'Position',[780,140,700,360]);
gridCL  = uigridlayout(panelCL,[1,2]);
axCorr  = uiaxes(gridCL); title(axCorr,'Correlation (window)');
axLag   = uiaxes(gridCL); title(axLag,'Lag (s, peak xcorr)');

% Optional: reduce UI right-click race conditions
% axVideo.ContextMenu = []; axCorr.ContextMenu = []; axLag.ContextMenu = [];
% for i = 1:numel(axTr), axTr(i).ContextMenu = []; end

% Bottom bar: slider row (1) + controls row (2)
bottom = uipanel(fig,'Position',[20,20,1460,100]);
g = uigridlayout(bottom,[2,16]);     % <-- was [2,14]
g.RowHeight = {36, 40};
g.ColumnWidth = {90,90,70,70,70,70,70, ...
                 '1x', ...           % slider stretch
                 90,70, ...         % (9-10)  Zoomed in + width
                 90,70, ...          % (11-12) Corr window label+spin 
                 90,70, ...         % (13-14) Max lag label+spin 
                 100,120};           % (15-16) auto + compute buttons
g.Padding = [8 6 8 6]; g.RowSpacing = 6; g.ColumnSpacing = 6;

% Row 1: time slider + label (unchanged)
timeSlider  = uislider(g,'ValueChangedFcn', @(s,~) updateTime(s.Value));
timeSlider.Layout.Row = 1; timeSlider.Layout.Column = [1 10];
timeLbl     = uilabel(g,'Text','Time: 0.00 s','HorizontalAlignment','right');
timeLbl.Layout.Row = 1; timeLbl.Layout.Column = 11;

% --- Refresh frequency control (Row 1) ---
lblRefresh = uilabel(g,'Text','Refresh freq (s):','HorizontalAlignment','right');
lblRefresh.Layout.Row    = 1; lblRefresh.Layout.Column = 13;

refreshSpin = uispinner(g,'Limits',[0.05 2],'Value',0.25,'Step',0.05, ...
    'Tooltip','Corr/Lag refresh interval (seconds)');
refreshSpin.Layout.Row    = 1; refreshSpin.Layout.Column = 14;


% Row 2: controls (left part unchanged)
csvBtn   = uibutton(g,'Text','Load CSV','ButtonPushedFcn', @(~,~) loadCSV());  csvBtn.Layout.Row=2; csvBtn.Layout.Column=1;
aviBtn   = uibutton(g,'Text','Load Video','Enable','off','ButtonPushedFcn', @(~,~) loadAVI());     aviBtn.Layout.Row=2; aviBtn.Layout.Column=2;
playBtn  = uibutton(g,'Text','Play','ButtonPushedFcn', @(~,~) playVideo());    playBtn.Layout.Row=2; playBtn.Layout.Column=3;
pauseBtn = uibutton(g,'Text','Pause','ButtonPushedFcn', @(~,~) pauseVideo());  pauseBtn.Layout.Row=2; pauseBtn.Layout.Column=4;
backBtn  = uibutton(g,'Text','<< 5s','ButtonPushedFcn', @(~,~) moveTime(-5));  backBtn.Layout.Row=2; backBtn.Layout.Column=5;
fwBtn    = uibutton(g,'Text','5s >>','ButtonPushedFcn', @(~,~) moveTime(5));   fwBtn.Layout.Row=2; fwBtn.Layout.Column=6;

smoothCB = uicheckbox(g,'Text','Smooth','Enable','off','ValueChangedFcn', @(~,~) onProcToggles());
smoothCB.Layout.Row=2; smoothCB.Layout.Column=7;

normCB   = uicheckbox(g,'Text','Normalize (ΔF/F%)','Enable','off','ValueChangedFcn', @(~,~) onProcToggles());
normCB.Layout.Row=2; normCB.Layout.Column=8;

% NEW: Zoom controls (go in the red-rectangle area)
zoomCB   = uicheckbox(g,'Text','Zoomed in','Enable','off', ...
    'Tooltip','Show only a time window centered at the current time', ...
    'ValueChangedFcn', @(~,~) applyZoom(timeSlider.Value));
zoomCB.Layout.Row=2; zoomCB.Layout.Column=9;

zoomSpin = uispinner(g,'Limits',[1 600],'Value',100,'Step',1,'Enable','off', ...
    'Tooltip','Zoom window width in seconds', ...
    'ValueChangedFcn', @(~,~) applyZoom(timeSlider.Value));
zoomSpin.Layout.Row=2; zoomSpin.Layout.Column=10;

% Correlation window (moved right)
lblWin   = uilabel(g,'Text','Corr window (s):','HorizontalAlignment','right');
lblWin.Layout.Row=2; lblWin.Layout.Column=11;
winSpin  = uispinner(g,'Limits',[0.5 120],'Value',60,'Step',0.5,'Enable','off', ...
    'ValueChangedFcn', @(~,~) requestCorrLagUpdate());
winSpin.Layout.Row=2; winSpin.Layout.Column=12;

% Max lag (moved right)
lblLag   = uilabel(g,'Text','Max lag (s):','HorizontalAlignment','right');
lblLag.Layout.Row=2; lblLag.Layout.Column=13;
lagSpin  = uispinner(g,'Limits',[0.05 60],'Value',30,'Step',0.05,'Enable','off', ...
    'ValueChangedFcn', @(~,~) requestCorrLagUpdate());
lagSpin.Layout.Row=2; lagSpin.Layout.Column=14;

autoCB   = uicheckbox(g,'Text','Auto update while playing','Enable','off');
autoCB.Layout.Row=2; autoCB.Layout.Column=15;

analyzeBtn = uibutton(g,'Text','Corr/Lag @ Now','Enable','off', ...
    'ButtonPushedFcn', @(~,~) computeCorrLagWindow(timeSlider.Value));
analyzeBtn.Layout.Row=2; analyzeBtn.Layout.Column=16;


%% Data/state
timeVec  = [];          % seconds
RawSig   = [];          % raw 470/410
roiNames = {};
aviObj   = [];
interpT  = [];          % frame-aligned times
vidTimer = timer('ExecutionMode','fixedRate','Period',1/30,'TimerFcn',@(~,~) playFrame());
nROI     = 0; fs = [];
lastCorrWall = tic;     % wall-clock throttle

% plot handles
hTrace = gobjects(1,8);         % trace lines
hVLine = gobjects(1,8);         % vertical cursor lines
yLims  = nan(8,2);              % cached y-limits
hVidIm = [];                    % video image handle

% preprocessing cache
procCache.data = [];
procCache.opt  = struct('smooth',false,'norm',false);
procCache.valid = false;

% --- circular connectivity window (created once, then reused)
circFig = [];      % classic figure (not uifigure) to avoid web UI races
circTL  = [];      % tiledlayout
circAxCorr = [];   % axes for correlation circular
circAxLag  = [];   % axes for lag/flow circular

%% --- Callbacks & helpers ------------------------------------------------
    function loadCSV()
        [f,p] = uigetfile({'*.csv','CSV files'},'Select CSV'); if isequal(f,0), return; end
        T = readtable(fullfile(p,f)); vars = T.Properties.VariableNames;

        ts = T{:,1};
        if isdatetime(ts), t = seconds(ts - ts(1));
        elseif isduration(ts), t = seconds(ts) - seconds(ts(1));
        else, t = ts - ts(1);
        end
        timeVec = t(:);
        fs = 1/mean(diff(timeVec));

        idx410 = find(contains(vars,'_sig_410','IgnoreCase',true));
        idx470 = find(contains(vars,'_sig_470','IgnoreCase',true));
        if isempty(idx410) || isempty(idx470)
            uialert(fig,'Expected "*_sig_410" and "*_sig_470" columns.','Error'); return;
        end
        nROI = min(8, min(numel(idx410), numel(idx470)));
        RawSig = zeros(length(timeVec),nROI);
        roiNames = cell(1,nROI);
        for k = 1:nROI
            RawSig(:,k) = T{:,idx470(k)} ./ T{:,idx410(k)};
            roiNames{k} = strrep(vars{idx410(k)},'_sig_410','');
        end

        % enable UI elements
        controls = [aviBtn, smoothCB, normCB, analyzeBtn, winSpin, lagSpin, autoCB, zoomCB, zoomSpin];
        arrayfun(@(h) set(h,'Enable','on'), controls);


        timeSlider.Limits = [0 max(timeVec)];
        timeSlider.Value  = 0; timeLbl.Text = 'Time: 0.00 s';

        procCache.valid = false;   % force preprocessing refresh
        initOrRefreshPlots();      % build persistent handles
        applyZoom(timeSlider.Value); 
        drawnow limitrate;
        cla(axCorr); cla(axLag);
    end

    function loadAVI()
        if isempty(timeVec), uialert(fig,'Load CSV first','Error'); return; end
        [f,p] = uigetfile({'*.avi;*.mp4','Video files'},'Select Video'); if isequal(f,0), return; end

        aviObj = VideoReader(fullfile(p,f));
        interpT = (0:1/aviObj.FrameRate:max(timeVec))';
        aviObj.CurrentTime = 0;

        if isvalid(vidTimer) && strcmp(vidTimer.Running,'on'), stop(vidTimer); end
        if isvalid(vidTimer), vidTimer.Period = 1/aviObj.FrameRate; end

        if hasFrame(aviObj)
            frame = readFrame(aviObj);
            if isempty(hVidIm) || ~isvalid(hVidIm)
                hVidIm = image(axVideo,'CData',frame); axis(axVideo,'off');
            else
                set(hVidIm,'CData',frame);
            end
        end
        drawnow limitrate;
    end

    function playVideo()
        if isempty(aviObj) || ~isvalid(vidTimer), return; end
        if ~strcmp(vidTimer.Running,'on'), start(vidTimer); end
    end

    function pauseVideo()
        if isvalid(vidTimer) && strcmp(vidTimer.Running,'on'), stop(vidTimer); end
    end

    function playFrame()
        % Defensive: stop if figure is gone
        if ~isvalid(fig), return; end

        % End-of-video or missing object
        if isempty(aviObj) || ~hasFrame(aviObj)
            if isvalid(vidTimer) && strcmp(vidTimer.Running,'on'), stop(vidTimer); end
            return;
        end

        % Update video image safely
        frame = readFrame(aviObj);
        if ~isempty(hVidIm) && isvalid(hVidIm) && isvalid(axVideo)
            set(hVidIm,'CData',frame);
        elseif isvalid(axVideo)
            hVidIm = image(axVideo,'CData',frame); axis(axVideo,'off');
        else
            return; % parent axis gone
        end

        % Advance time / cursor
        tnow = aviObj.CurrentTime;
        if isvalid(fig), updateTime(tnow, true); end

        % Throttle corr/lag by wall-clock (by user setting)
        if isvalid(fig) && autoCB.Value && toc(lastCorrWall) >= refreshSpin.Value
            if isvalid(axCorr) && isvalid(axLag)
                computeCorrLagWindow(tnow);
            end
            lastCorrWall = tic;
        end
        drawnow limitrate;
    end

    function moveTime(delta)
        t0 = min(max(timeSlider.Value + delta, timeSlider.Limits(1)), timeSlider.Limits(2));
        updateTime(t0, false);
        if autoCB.Value, requestCorrLagUpdate(); end
    end

    function updateTime(t, varargin)
        fromTimer = false; if ~isempty(varargin), fromTimer = logical(varargin{1}); end
        if ~isvalid(timeSlider) || ~isvalid(timeLbl), return; end
        timeSlider.Value = t; timeLbl.Text = sprintf('Time: %.2f s',t);

        % update vertical cursor only
        for m = 1:nROI
            if isgraphics(hVLine(m)), set(hVLine(m),'XData',[t t]); end
        end

        % NEW: update zoomed view if enabled
        applyZoom(t);

        % Immediate recompute on user scrub
        if ~fromTimer && autoCB.Value
            computeCorrLagWindow(t);
            lastCorrWall = tic;   % reset throttle
        end

        % When scrubbing, seek video to nearest frame and show it
        if ~fromTimer && ~isempty(aviObj) && ~isempty(interpT)
            idx = find(interpT >= t,1);
            if ~isempty(idx)
                aviObj.CurrentTime = interpT(idx);
                if hasFrame(aviObj)
                    frame = readFrame(aviObj);
                    if ~isempty(hVidIm) && isvalid(hVidIm)
                        set(hVidIm,'CData',frame);
                    else
                        hVidIm = image(axVideo,'CData',frame); axis(axVideo,'off');
                    end
                end
            end
        end
    end

    function onProcToggles()
        procCache.valid = false;
        refreshTraceYData();
        requestCorrLagUpdate();
    end

    function data = getProcessedData()
        if procCache.valid ...
                && procCache.opt.smooth == logical(smoothCB.Value) ...
                && procCache.opt.norm   == logical(normCB.Value)
            data = procCache.data; return;
        end
        data = RawSig; if isempty(data), return; end

        if smoothCB.Value
            data = smoothdata(data,'movmean',5);
        end
        if normCB.Value
            % robust baseline: 10th percentile per ROI
            F0 = prctile(data,10,1);
            F0(F0==0) = eps;
            data = bsxfun(@rdivide, bsxfun(@minus, data, F0), F0) * 100; % ΔF/F%
        end

        procCache.data = data;
        procCache.opt  = struct('smooth',logical(smoothCB.Value),'norm',logical(normCB.Value));
        procCache.valid = true;
    end

    function initOrRefreshPlots()
        data = getProcessedData();
        for m = 1:8
            cla(axTr(m)); hold(axTr(m),'on');
            if m <= nROI
                hTrace(m) = plot(axTr(m), timeVec, data(:,m), 'LineWidth',1);
                yl = [min(data(:,m)) max(data(:,m))];
                if yl(1)==yl(2), yl = yl + [-1 1]*0.5; end
                yLims(m,:) = yl; ylim(axTr(m),yl);
                hVLine(m) = plot(axTr(m), [timeSlider.Value timeSlider.Value], yl, 'k--');
                title(axTr(m), roiNames{m});
                if normCB.Value, ylabel(axTr(m),'ΔF/F (%)'); else, ylabel(axTr(m),'470/410'); end
            else
                title(axTr(m),''); hTrace(m) = gobjects(1); hVLine(m) = gobjects(1);
            end
            hold(axTr(m),'off');
            if m < 8, axTr(m).XTick = []; axTr(m).XColor='none'; end
        end
        xlabel(axTr(8),'Time (s)'); axTr(8).XColor='k'; axTr(8).XTickMode='auto';
    end

    function refreshTraceYData()
        if isempty(RawSig) || any(~isgraphics(hTrace(1:nROI))), return; end
        data = getProcessedData();
        for m = 1:nROI
            set(hTrace(m),'YData',data(:,m));
            % keep stable y-lims unless change >10%
            cur = yLims(m,:);
            new = [min(data(:,m)) max(data(:,m))];
            if new(1)==new(2), new = new + [-1 1]*0.5; end
            if any(~isfinite(cur)) || max(abs((new-cur)./max(1,abs(cur)))) > 0.1
                yLims(m,:) = new; ylim(axTr(m),new);
                if isgraphics(hVLine(m)), set(hVLine(m),'YData',new); end
            end
            if normCB.Value, ylabel(axTr(m),'ΔF/F (%)'); else, ylabel(axTr(m),'470/410'); end
        end
        drawnow limitrate;
    end

    function applyZoom(tCenter)
        if isempty(timeVec) || ~isvalid(zoomCB) || ~zoomCB.Value
            % Show full time range
            for m = 1:nROI
                if isgraphics(axTr(m))
                    xlim(axTr(m), [timeSlider.Limits(1) timeSlider.Limits(2)]);
                    ylim(axTr(m), yLims(m,:));   % restore cached full-range limits
                end
            end
            return;
        end

        W  = max(1, zoomSpin.Value);  % total width (s)
        hw = W/2;
        t1 = max(timeSlider.Limits(1), tCenter - hw);
        t2 = min(timeSlider.Limits(2), tCenter + hw);

        % Expand window if near edges
        if (t2 - t1) < W
            if t1 <= timeSlider.Limits(1)
                t2 = min(timeSlider.Limits(2), t1 + W);
            elseif t2 >= timeSlider.Limits(2)
                t1 = max(timeSlider.Limits(1), t2 - W);
            end
        end

        data = getProcessedData();

        for m = 1:nROI
            if isgraphics(axTr(m))
                % set X window
                xlim(axTr(m), [t1 t2]);

                % NEW: set Y limits only to min/max within [t1,t2]
                idx = (timeVec >= t1 & timeVec <= t2);
                if any(idx)
                    ymin = min(data(idx,m));
                    ymax = max(data(idx,m));
                    if ymin==ymax, ymin=ymin-0.5; ymax=ymax+0.5; end
                    ylim(axTr(m), [ymin ymax]);
                    if isgraphics(hVLine(m)), set(hVLine(m),'YData',[ymin ymax]); end
                end
            end
        end
    end




    function requestCorrLagUpdate()
        if isempty(RawSig), return; end
        if autoCB.Value
            % timer will handle throttled updates
            return;
        else
            computeCorrLagWindow(timeSlider.Value);
        end
    end

    function computeCorrLagWindow(tCenter)
        if isempty(RawSig), uialert(fig,'Load data first','Error'); return; end
        if ~isvalid(fig) || ~isvalid(axCorr) || ~isvalid(axLag), return; end

        win = winSpin.Value;          % seconds (half-width)
        maxLagSec = lagSpin.Value;    % seconds

        t1 = max(timeSlider.Limits(1), tCenter - win);
        t2 = min(timeSlider.Limits(2), tCenter + win);

        idx = (timeVec >= t1) & (timeVec <= t2);
        if nnz(idx) < 5
            cla(axCorr); cla(axLag);
            title(axCorr,sprintf('Correlation (%.2f-%.2f s): insufficient samples',t1,t2));
            title(axLag,''); return;
        end

        dataW = getProcessedData();
        dataW = dataW(idx,1:nROI);

        % correlation
        % correlation
        C = corrcoef(dataW);

        % --- Lag matrix with proper caps and |xcorr| peak ---
        % Optional downsample for speed:
        % --- Lag matrix with proper caps and |xcorr| peak (robust near-zero tie-break) ---
        ds = 1;                              % keep at 1 for your 20 Hz data
        if ds > 1, dataW = resample(dataW,1,ds); end
        fs_eff = fs / ds;

        dataZ  = dataW - mean(dataW,1);
        Nsamp  = size(dataZ,1);
        maxLagSec_req = lagSpin.Value;
        maxLag       = max(1, min(round(maxLagSec_req*fs_eff), Nsamp-1));

        lagM = zeros(nROI);
        for i = 1:nROI
            xi = dataZ(:,i);
            for j = i+1:nROI
                xj = dataZ(:,j);

                % Prefer finddelay if available (positive => xj lags xi)
                try
                    D = finddelay(xi, xj);           % Signal Processing Toolbox
                    if abs(D) > maxLag               % cap by requested window
                        D = sign(D) * maxLag;
                    end
                    lagSec = D / fs_eff;
                catch
                    % Cross-corr with near-peak tie-break toward zero lag
                    [xc,lags] = xcorr(xi,xj,maxLag,'coeff');
                    amax = max(abs(xc));
                    tol  = 0.02;                     % accept within 2% of the max
                    cand = find(abs(xc) >= (1-tol)*amax);

                    % among near-peak candidates, choose the lag with smallest |lag|
                    [~,kmin] = min(abs(lags(cand)));
                    ix = cand(kmin);

                    lagSec = lags(ix) / fs_eff;      % + => j lags i
                end

                lagM(i,j) = lagSec;
                lagM(j,i) = -lagSec;
            end
        end

        % draw
        imagesc(axCorr,C,[-1 1]); axis(axCorr,'square'); colorbar(axCorr);
        xticks(axCorr,1:nROI); yticks(axCorr,1:nROI);
        xticklabels(axCorr,roiNames(1:nROI)); yticklabels(axCorr,roiNames(1:nROI));
        title(axCorr,sprintf('Correlation (%.2f to %.2f s)',t1,t2));

        imagesc(axLag,lagM); axis(axLag,'square'); colorbar(axLag);
        xticks(axLag,1:nROI); yticks(axLag,1:nROI);
        xticklabels(axLag,roiNames(1:nROI)); yticklabels(axLag,roiNames(1:nROI));
        effLagSec = maxLag/fs_eff;
        lagNote = '';
        if effLagSec < maxLagSec_req - 1e-9
            lagNote = sprintf(' (capped at %.2fs by window)', effLagSec);
        end
        title(axLag,sprintf('Lag (s), max |xcorr| (maxLag=%.2fs)%s', maxLagSec_req, lagNote));
        drawnow limitrate;

        updateCircular(C, lagM);



    end

%% Keyboard shortcuts
fig.KeyPressFcn = @(~,evt) onKey(evt);
    function onKey(evt)
        switch lower(evt.Key)
            case 'c', computeCorrLagWindow(timeSlider.Value);
            case 'space'
                if isvalid(vidTimer) && strcmp(vidTimer.Running,'on'), pauseVideo(); else, playVideo(); end
        end
    end
    function ensureCircularWindow()
        if isempty(circFig) || ~isvalid(circFig)
            circFig = figure('Name','Circular connectivity', ...
                'Color','w','Position',[100 100 780 360]);
            circTL = tiledlayout(circFig,1,2,'Padding','compact','TileSpacing','compact');
            circAxCorr = nexttile(circTL,1); title(circAxCorr,'Correlation (circular)');
            circAxLag  = nexttile(circTL,2); title(circAxLag,'Lag/flow (circular)');
        end
    end

    function updateCircular(C, lagM)
        % draw into the existing axes (no new figures)
        ensureCircularWindow();
        if isempty(roiNames), return; end
        names = roiNames(1:nROI);

        % ---- Correlation (undirected, continuous colors) ----
        dens = 0.35;                               % keep top 35% |r|
        A = C; A(1:nROI+1:end) = 0;
        W = abs(A);
        cla(circAxCorr);

        if any(W(:))
            wvec = W(triu(true(nROI),1));
            thr  = quantile(wvec(wvec>0), max(0,1-dens));
            mask = triu(W>=thr,1);
            [iu,ju] = find(mask);
            N  = nROI;

            % Build graph (keep all N nodes so labels length matches)
            Gu = graph(iu, ju, [], N);
            h1 = plot(circAxCorr, Gu,'Layout','circle','NodeLabel',names, ...
                'MarkerSize',8,'NodeColor',[.15 .15 .15],'LineWidth',1.5);

            if ~isempty(iu)
                % Edge stats
                wu   = A(sub2ind([nROI nROI],iu,ju));  % signed r for drawn edges
                wabs = abs(wu);

                % Line width by |r|
                lw = 0.5 + 4*(wabs - min(wabs))/(max(wabs)-min(wabs)+eps);
                h1.LineWidth = lw;

                % --- Continuous edge colors by r (blue→white→red) ---
                % robust symmetric scale around 0, so white ≈ 0
                cmax = max(0.4, max(abs(wu)));     % clamp to avoid over-saturation
                cmapSize = 256;
                % build a simple diverging colormap [blue; white; red]
                cmapSize = 256;
                stops = [0 0.5 1];
                cols  = [0 0.2 0.8; 1 1 1; 0.9 0.1 0.1];   % blue—white—red
                xi    = linspace(0,1,cmapSize)';
                cmapD = interp1(stops, cols, xi);

                colormap(circAxCorr, cmapD);
                caxis(circAxCorr, [-1 1]);                 % <<< fixed range for r
                h1.EdgeCData = wu;                          % color edges by r
                cb = colorbar(circAxCorr); cb.Label.String = 'Correlation r';

                % % If your MATLAB supports EdgeCData (R2018b+), use true colormapping:
                % if isprop(h1,'EdgeCData')
                %     h1.EdgeCData = wu;                     % one value per edge
                %     % colormap(circAxCorr, cmapD);
                %     colormap(circAxCorr, [-1 1]);
                %     % caxis(circAxCorr, [-cmax cmax]);       % symmetric limits
                %     cb = colorbar(circAxCorr);
                %     cb.Label.String = 'Correlation r';
                % else
                %     % Fallback: map r → RGB manually, still show a legend colorbar
                %     t   = (wu + cmax) ./ (2*cmax + eps);   % map r∈[-cmax,cmax] → [0,1]
                %     idx = max(1, min(cmapSize, round(t*(cmapSize-1))+1));
                %     h1.EdgeColor = cmapD(idx,:);
                %     colormap(circAxCorr, cmapD);
                %     caxis(circAxCorr, [-cmax cmax]);
                %     cb = colorbar(circAxCorr); cb.Label.String = 'Correlation r';
                % end
            else
                axis(circAxCorr,'off'); text(circAxCorr,0.5,0.5,'No edges','Horiz','center');
            end
        end

        title(circAxCorr, sprintf('Correlation circular (top %.0f%% |r|)', dens*100));


        % ---- Lag/flow (directed) ----
        % ---- Lag/flow (directed) ----
lagThr = 0.3;                           % show 0.5, 0.8, 1.2, 2.2 s delays
cla(circAxLag);

leadMask = abs(lagM) >= lagThr;
[iU,jU]  = find(triu(leadMask,1));      % use upper triangle once
src=[]; dst=[]; lagv=[]; rv=[];

for k=1:numel(iU)
    i=iU(k); j=jU(k);
    if lagM(i,j) > 0      % j lags i  => edge i -> j
        src(end+1)=i; dst(end+1)=j; lagv(end+1)=lagM(i,j); rv(end+1)=C(i,j);
    elseif lagM(i,j) < 0  % i lags j  => edge j -> i
        src(end+1)=j; dst(end+1)=i; lagv(end+1)=abs(lagM(i,j)); rv(end+1)=C(i,j);
    end
end

if isempty(src)
    axis(circAxLag,'off'); title(circAxLag,'Lag/flow (no edges ≥ threshold)');
else
    Gd = digraph(src, dst, [], nROI);
    h2 = plot(circAxLag, Gd,'Layout','circle','NodeLabel',names, ...
              'MarkerSize',8,'NodeColor',[.15 .15 .15],'ArrowSize',12);
    % edge width by |r|
    rabs = max(0, min(1, abs(rv)));
    w2 = 0.5 + 4*(rabs - min(rabs))/(max(rabs)-min(rabs)+eps);
    h2.LineWidth = w2;

    % edge color by lag magnitude
    L  = lagv(:);
    Ln = (L - min(L)) / (max(L)-min(L) + eps);
    cmapLag = parula(256);
    idx = max(1, round(Ln*255)+1);
    h2.EdgeColor = cmapLag(idx,:);
    colormap(circAxLag, cmapLag);
    cb = colorbar(circAxLag); cb.Label.String = 'Lag (s)  leader → follower';
    title(circAxLag, sprintf('Lag/flow (|lag| ≥ %.1fs)', lagThr));

    % OPTIONAL: print lag value on each edge for sanity‑check
    try
        for e = 1:numedges(Gd)
            mid = mean([h2.XData(src(e)) h2.XData(dst(e)); h2.YData(src(e)) h2.YData(dst(e))],2);
            text(circAxLag, mid(1), mid(2), sprintf('%.1fs', lagv(e)), ...
                 'HorizontalAlignment','center','FontSize',8,'Color',[0 0 0]);
        end
    end
end

    end

%% ---- Cleanup -----------------------------------------------------------
    function cleanupAndClose()
        isClosing = true; drawnow;  % let pending UI events settle
        try
            if exist('vidTimer','var') && isvalid(vidTimer)
                stop(vidTimer);
                set(vidTimer,'TimerFcn',[]);   % detach callback
                delete(vidTimer);
            end
        end
        try
            if exist('fig','var') && isvalid(fig)
                % Detach all app callbacks to avoid late hits
                fig.KeyPressFcn = [];
                timeSlider.ValueChangedFcn = [];
                % (do same for any other UI callbacks you created)
                delete(fig);
            end
        end
        try
            if ~isempty(circFig) && isvalid(circFig), delete(circFig); end
        end
    end

end
