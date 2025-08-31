function fiber_photometry_behavior_analysis_combined_sorted_V1_3()
% Fiber Photometry: Multi-file analysis with Aggregate 2x2 and Per-mouse 4x4 pages
% - Imports multiple CSV (raw) + Excel (Behavior1/Behavior2 tags) pairs
% - Smooths Real_Signal (470/410) with LP Butterworth
% - Aligns to tags, computes ΔF/F, then scales to PERCENT (%)
% - Aggregates across mice and sorts trials by AUC in user-defined windows
% - Final aggregate figure: 2x2 (B1 overlay, B1 heatmap, B2 overlay, B2 heatmap)
% - Per-mouse pages: 4x4, each row = one mouse: B1 traces, B1 HM, B2 traces, B2 HM
%
% Version: V1.3
% Change log:
%   - ΔF/F now plotted as percent everywhere (traces & heatmaps)
%   - Heatmaps sorted by AUC in ASCENDING order (lowest→highest)
%   - Matching sort used for overlays by default (toggle via sortOverlayLikeHeatmap)

%% ===== Settings =====
close all; clearvars; clc; warning off;

% Display/Global limits — now in PERCENT units
Global_heatmap_Lim = [0 10];     % e.g., 0–20%
Global_Y_axis_Lim  = [-2 5];     % e.g., -2%–+5%

% Analysis window around onset (sec)
preTime  = -5;     % seconds before tag
postTime =  5;     % seconds after tag

% AUC windows (sec) for sorting
aucWindow_B1 = [0,  min(100, postTime)];
aucWindow_B2 = [-1, 3];

% Sorting orders
sortOrderHeatmap = 'ascend';  % user request: heatmaps ascend (lowest at top)
sortOverlayLikeHeatmap = true; % if true, overlays use same sort order for consistency

% Figure shapes
makePerMouseFigureSquare = true;  % square overall window for per-mouse pages

%% ===== File selection =====
[csvFileNames, csvPath] = uigetfile('*.csv', 'Select the Raw CSV Data Files', 'MultiSelect', 'on');
if isequal(csvFileNames,0), disp('User canceled CSV selection.'); return; end
if ischar(csvFileNames), csvFileNames = {csvFileNames}; end
numCsvFiles = numel(csvFileNames);

[xlsxFileNames, xlsxPath] = uigetfile('*.xlsx', 'Select the Behavior Tags Excel Files', 'MultiSelect', 'on');
if isequal(xlsxFileNames,0), disp('User canceled Excel selection.'); return; end
if ischar(xlsxFileNames), xlsxFileNames = {xlsxFileNames}; end
numXlsxFiles = numel(xlsxFileNames);

if numCsvFiles ~= numXlsxFiles
    error('Number of CSV files (%d) must match number of Excel files (%d).', numCsvFiles, numXlsxFiles);
end
fprintf('Number of pairs: %d\n', numCsvFiles);

%% ===== Aggregation containers =====
all_normalizedData1 = []; all_timeWindow1 = [];
all_normalizedData2 = []; all_timeWindow2 = [];

% Per-mouse containers
perMouse = struct('name',{},'init',{}, ...
                  't1',{},'norm1',{},'mean1',{},'sem1',{}, ...
                  't2',{},'norm2',{},'mean2',{},'sem2',{});

% Helper for initials from filename
get_initials = @(fname) strtrim(regexprep( ...
                      upper(regexprep(regexprep(fname,'\.[^.]*$',''), ...
                      '(^|[^A-Za-z0-9])([A-Za-z0-9])[A-Za-z0-9]*',' $2')), ...
                      ' +',''));
trim_initials = @(s) s(1:min(8,numel(s)));

%% ===== Main loop: per file pair =====
for i = 1:numCsvFiles
    fprintf('\nProcessing pair %d/%d\n', i, numCsvFiles);
    csvFileName  = csvFileNames{i};
    xlsxFileName = xlsxFileNames{i};
    fprintf('CSV:  %s\n', csvFileName);
    fprintf('XLSX: %s\n', xlsxFileName);

    csvFullPath  = fullfile(csvPath,  csvFileName);
    xlsxFullPath = fullfile(xlsxPath, xlsxFileName);

    % -- Load CSV
    try
        rawData = readtable(csvFullPath);
    catch ME
        warning('Failed reading CSV %s: %s. Skipping.', csvFileName, ME.message);
        continue;
    end
    if width(rawData) < 3
        warning('CSV %s lacks required first 3 columns. Skipping.', csvFileName);
        continue;
    end

    % -- Extract and validate signal
    TimeStamp = rawData{:,1};
    LED_410   = rawData{:,2};
    LED_470   = rawData{:,3};

    if isduration(TimeStamp)
        timeVector = seconds(TimeStamp) - seconds(TimeStamp(1));
    elseif isdatetime(TimeStamp)
        timeVector = seconds(TimeStamp - TimeStamp(1));
    else
        warning('Timestamp in %s not datetime/duration. Skipping.', csvFileName);
        continue;
    end
    timeVector = timeVector(:);
    Real_Signal = (LED_470(:) ./ LED_410(:));

    finiteIdx = isfinite(timeVector) & isfinite(Real_Signal);
    if ~all(finiteIdx)
        timeVector = timeVector(finiteIdx);
        Real_Signal = Real_Signal(finiteIdx);
    end
    if isempty(timeVector) || isempty(Real_Signal)
        warning('No valid data in %s after cleaning. Skipping.', csvFileName);
        continue;
    end

    % -- Smooth with LP Butterworth
    Fc = 5; order = 4;
    Fs = 1 / median(diff(timeVector));
    [b, a] = butter(order, Fc/(Fs/2), 'low');
    Real_Signal_Smoothed = filtfilt(b, a, Real_Signal);

    % -- Load behavior tags
    behavior1Times = []; behavior2Times = [];
    try
        [~, sheets] = xlsfinfo(xlsxFullPath);
        if ismember('Behavior1', sheets)
            T1 = readtable(xlsxFullPath, 'Sheet','Behavior1');
            if ismember('Behavior1_Onset_Time_s', T1.Properties.VariableNames)
                behavior1Times = T1.Behavior1_Onset_Time_s;
            end
        end
        if ismember('Behavior2', sheets)
            T2 = readtable(xlsxFullPath, 'Sheet','Behavior2');
            if ismember('Behavior2_Onset_Time_s', T2.Properties.VariableNames)
                behavior2Times = T2.Behavior2_Onset_Time_s;
            end
        end
    catch ME
        warning('Failed reading tags %s: %s', xlsxFileName, ME.message);
    end

    % -- Compute aligned and ΔF/F for each behavior (ΔF/F returns PERCENT)
    thisMouse.name = csvFileName;
    thisMouse.init = trim_initials(get_initials(csvFileName));

    % Behavior 1
    if ~isempty(behavior1Times)
        [aligned1, t1] = extract_aligned_traces(timeVector, Real_Signal_Smoothed, behavior1Times, preTime, postTime);
        if ~isempty(aligned1)
            norm1 = compute_deltaF_over_F_percent(aligned1, t1, preTime, postTime); % <-- percent
            [m1, s1] = compute_mean_sem(norm1);
            % aggregate
            all_normalizedData1 = [all_normalizedData1; norm1];
            if isempty(all_timeWindow1), all_timeWindow1 = t1; end
            % per-mouse
            thisMouse.t1 = t1(:)'; thisMouse.norm1 = norm1; thisMouse.mean1 = m1; thisMouse.sem1 = s1;
        else
            thisMouse.t1=[]; thisMouse.norm1=[]; thisMouse.mean1=[]; thisMouse.sem1=[];
        end
    else
        thisMouse.t1=[]; thisMouse.norm1=[]; thisMouse.mean1=[]; thisMouse.sem1=[];
    end

    % Behavior 2
    if ~isempty(behavior2Times)
        [aligned2, t2] = extract_aligned_traces(timeVector, Real_Signal_Smoothed, behavior2Times, preTime, postTime);
        if ~isempty(aligned2)
            norm2 = compute_deltaF_over_F_percent(aligned2, t2, preTime, postTime); % <-- percent
            [m2, s2] = compute_mean_sem(norm2);
            % aggregate
            all_normalizedData2 = [all_normalizedData2; norm2];
            if isempty(all_timeWindow2), all_timeWindow2 = t2; end
            % per-mouse
            thisMouse.t2 = t2(:)'; thisMouse.norm2 = norm2; thisMouse.mean2 = m2; thisMouse.sem2 = s2;
        else
            thisMouse.t2=[]; thisMouse.norm2=[]; thisMouse.mean2=[]; thisMouse.sem2=[];
        end
    else
        thisMouse.t2=[]; thisMouse.norm2=[]; thisMouse.mean2=[]; thisMouse.sem2=[];
    end

    perMouse(end+1) = thisMouse; %#ok<SAGROW>
end

%% ===== Aggregate sorting and 2x2 figure =====
% Behavior 1 sorting (by AUC within window)
sorted1 = [];
if ~isempty(all_normalizedData1)
    idxAUC1 = all_timeWindow1 >= aucWindow_B1(1) & all_timeWindow1 <= aucWindow_B1(2);
    if any(idxAUC1)
        aucVals1 = trapz(all_timeWindow1(idxAUC1), all_normalizedData1(:,idxAUC1), 2);
        [~, ord1] = sort(aucVals1, sortOrderHeatmap); % ASCEND for heatmap
    else
        ord1 = (1:size(all_normalizedData1,1)).';
    end
    sorted1 = all_normalizedData1(ord1,:);
    mean1   = mean(sorted1,1,'omitnan');
    sem1    = std(sorted1,0,1,'omitnan')/sqrt(size(sorted1,1));
end

% Behavior 2 sorting
sorted2 = [];
if ~isempty(all_normalizedData2)
    idxAUC2 = all_timeWindow2 >= aucWindow_B2(1) & all_timeWindow2 <= aucWindow_B2(2);
    if any(idxAUC2)
        aucVals2 = trapz(all_timeWindow2(idxAUC2), all_normalizedData2(:,idxAUC2), 2);
        [~, ord2] = sort(aucVals2, sortOrderHeatmap); % ASCEND for heatmap
    else
        ord2 = (1:size(all_normalizedData2,1)).';
    end
    sorted2 = all_normalizedData2(ord2,:);
    mean2   = mean(sorted2,1,'omitnan');
    sem2    = std(sorted2,0,1,'omitnan')/sqrt(size(sorted2,1));
end

% ---- Aggregate 2x2 ----
figAgg = figure('Name','Aggregate Behaviors (2x2)','NumberTitle','off','Position',[10 10 1100 900]);

% TL: B1 overlay (stack of trials + mean±SEM)
subplot(2,2,1);
if ~isempty(sorted1)
    hold on;
    toPlot1 = sorted1; % (sorted to match heatmap)
    plot(all_timeWindow1, toPlot1', 'Color',[0.8 0.8 0.8], 'LineWidth',1);
    plot(all_timeWindow1, mean1, 'k','LineWidth',2);
    x_fill = [all_timeWindow1(:); flipud(all_timeWindow1(:))];
    y_fill = [mean1(:)+sem1(:); flipud(mean1(:)-sem1(:))];
    fill(x_fill, y_fill, 'k','FaceAlpha',0.3,'EdgeColor','none');
    xlabel('Time (s)'); ylabel('\DeltaF/F (%)'); title('Aggregate B1 (AUC-sorted ↑)');
    grid on; ylim(Global_Y_axis_Lim); box off;
    axis('square');
else
    text(0.5,0.5,'No B1 data','HorizontalAlignment','center'); axis off;
end

% TR: B1 heatmap
subplot(2,2,2);
if ~isempty(sorted1)
    imagesc(all_timeWindow1, 1:size(sorted1,1), sorted1);
    set(gca,'YDir','normal'); colormap('parula'); caxis(Global_heatmap_Lim);
    cb = colorbar; ylabel(cb,'\DeltaF/F (%)');
    xlabel('Time (s)'); ylabel('Trial'); title('B1 Heatmap (AUC ↑)');
else
    text(0.5,0.5,'No B1 data','HorizontalAlignment','center'); axis off;
end

% BL: B2 overlay
subplot(2,2,3);
if ~isempty(sorted2)
    hold on;
    toPlot2 = sorted2; % (sorted to match heatmap)
    plot(all_timeWindow2, toPlot2', 'Color',[0.8 0.8 0.8], 'LineWidth',1);
    plot(all_timeWindow2, mean2, 'k','LineWidth',2);
    x_fill = [all_timeWindow2(:); flipud(all_timeWindow2(:))];
    y_fill = [mean2(:)+sem2(:); flipud(mean2(:)-sem2(:))];
    fill(x_fill, y_fill, 'k','FaceAlpha',0.3,'EdgeColor','none');
    xlabel('Time (s)'); ylabel('\DeltaF/F (%)'); title('Aggregate B2 (AUC-sorted ↑)');
    grid on; ylim(Global_Y_axis_Lim); box off;
    axis('square');
else
    text(0.5,0.5,'No B2 data','HorizontalAlignment','center'); axis off;
end

% BR: B2 heatmap
subplot(2,2,4);
if ~isempty(sorted2)
    imagesc(all_timeWindow2, 1:size(sorted2,1), sorted2);
    set(gca,'YDir','normal'); colormap('parula'); caxis(Global_heatmap_Lim);
    cb = colorbar; ylabel(cb,'\DeltaF/F (%)');
    xlabel('Time (s)'); ylabel('Trial'); title('B2 Heatmap (AUC ↑)');
else
    text(0.5,0.5,'No B2 data','HorizontalAlignment','center'); axis off;
end

%% ===== Per-mouse 4x4 pages (4 mice per figure) =====
nMice = numel(perMouse);
if nMice == 0
    warning('No per-mouse data to plot.');
else
    micePerPage = 4;
    nPages = ceil(nMice / micePerPage);

    for pg = 1:nPages
        iStart = (pg-1)*micePerPage + 1;
        iEnd   = min(pg*micePerPage, nMice);

        if makePerMouseFigureSquare
            figPM = figure('Name',sprintf('Per-mouse Behaviors (page %d/%d)',pg,nPages), ...
                           'NumberTitle','off','Position',[20 20 1000 1000]); % square
        else
            figPM = figure('Name',sprintf('Per-mouse Behaviors (page %d/%d)',pg,nPages), ...
                           'NumberTitle','off','Position',[20 20 1300 900]');
        end

        for r = 1:(iEnd - iStart + 1)
            m = perMouse(iStart + r - 1);

            c1 = (r-1)*4 + 1;  % B1 traces
            c2 = (r-1)*4 + 2;  % B1 heatmap
            c3 = (r-1)*4 + 3;  % B2 traces
            c4 = (r-1)*4 + 4;  % B2 heatmap

            % --- B1 traces
            ax1 = subplot(4,4,c1); hold(ax1,'on');
            if ~isempty(m.t1) && ~isempty(m.norm1)
                % local sort for matching HM
                idxB1 = m.t1 >= aucWindow_B1(1) & m.t1 <= aucWindow_B1(2);
                if any(idxB1)
                    aucLocal1 = trapz(m.t1(idxB1), m.norm1(:,idxB1), 2);
                    [~, o1] = sort(aucLocal1, sortOrderHeatmap); % ascend for the HM
                else
                    o1 = (1:size(m.norm1,1)).';
                end
                if sortOverlayLikeHeatmap
                    dat1_traces = m.norm1(o1,:);
                else
                    dat1_traces = m.norm1;
                end

                plot(m.t1, dat1_traces', 'Color',[0.85 0.85 0.85], 'LineWidth',0.8);
                % recompute mean/sem for the plotted order (not necessary for mean curve)
                if isempty(m.mean1) || sortOverlayLikeHeatmap
                    mean1m = mean(dat1_traces,1,'omitnan');
                    sem1m  = std(dat1_traces,0,1,'omitnan')/sqrt(size(dat1_traces,1));
                else
                    mean1m = m.mean1; sem1m = m.sem1;
                end
                plot(m.t1, mean1m, 'k', 'LineWidth',1.5);
                x_fill = [m.t1(:); flipud(m.t1(:))];
                y_fill = [mean1m(:)+sem1m(:); flipud(mean1m(:)-sem1m(:))];
                fill(x_fill, y_fill, 'k', 'FaceAlpha',0.25,'EdgeColor','none');

                ylim(Global_Y_axis_Lim);
                title(sprintf('%s – B1', m.init),'FontWeight','normal');
            else
                text(0.5,0.5,'No B1','HorizontalAlignment','center'); axis(ax1,'off');
            end
            grid on; box off; xlabel(''); ylabel('\DeltaF/F (%)');
            shrink_axes(ax1);

            % --- B1 heatmap (AUC ascend)
            ax2 = subplot(4,4,c2);
            if ~isempty(m.t1) && ~isempty(m.norm1)
                if ~exist('o1','var') || isempty(o1), o1 = (1:size(m.norm1,1)).'; end
                dat1_hm = m.norm1(o1,:);
                imagesc(m.t1, 1:size(dat1_hm,1), dat1_hm); set(ax2,'YDir','normal');
                colormap(ax2,'parula'); caxis(ax2,Global_heatmap_Lim);
                cb = colorbar(ax2); ylabel(cb,'\DeltaF/F (%)');
                title('B1 HM (AUC ↑)','FontWeight','normal'); xlabel(''); ylabel('');
            else
                text(0.5,0.5,'No B1','HorizontalAlignment','center'); axis(ax2,'off');
            end
            axis('square');
            clear o1;

            % --- B2 traces
            ax3 = subplot(4,4,c3); hold(ax3,'on');
            if ~isempty(m.t2) && ~isempty(m.norm2)
                idxB2 = m.t2 >= aucWindow_B2(1) & m.t2 <= aucWindow_B2(2);
                if any(idxB2)
                    aucLocal2 = trapz(m.t2(idxB2), m.norm2(:,idxB2), 2);
                    [~, o2] = sort(aucLocal2, sortOrderHeatmap); % ascend for the HM
                else
                    o2 = (1:size(m.norm2,1)).';
                end
                if sortOverlayLikeHeatmap
                    dat2_traces = m.norm2(o2,:);
                else
                    dat2_traces = m.norm2;
                end

                plot(m.t2, dat2_traces', 'Color',[0.85 0.85 0.85], 'LineWidth',0.8);
                if isempty(m.mean2) || sortOverlayLikeHeatmap
                    mean2m = mean(dat2_traces,1,'omitnan');
                    sem2m  = std(dat2_traces,0,1,'omitnan')/sqrt(size(dat2_traces,1));
                else
                    mean2m = m.mean2; sem2m = m.sem2;
                end
                plot(m.t2, mean2m, 'k', 'LineWidth',1.5);
                x_fill = [m.t2(:); flipud(m.t2(:))];
                y_fill = [mean2m(:)+sem2m(:); flipud(mean2m(:)-sem2m(:))];
                fill(x_fill, y_fill, 'k', 'FaceAlpha',0.25,'EdgeColor','none');

                ylim(Global_Y_axis_Lim);
                title(sprintf('%s – B2', m.init),'FontWeight','normal');
            else
                text(0.5,0.5,'No B2','HorizontalAlignment','center'); axis(ax3,'off');
            end
            grid on; box off; xlabel(''); ylabel('\DeltaF/F (%)');
            shrink_axes(ax3);

            % --- B2 heatmap (AUC ascend)
            ax4 = subplot(4,4,c4);
            if ~isempty(m.t2) && ~isempty(m.norm2)
                if ~exist('o2','var') || isempty(o2), o2 = (1:size(m.norm2,1)).'; end
                dat2_hm = m.norm2(o2,:);
                imagesc(m.t2, 1:size(dat2_hm,1), dat2_hm); set(ax4,'YDir','normal');
                colormap(ax4,'parula'); caxis(ax4,Global_heatmap_Lim);
                cb = colorbar(ax4); ylabel(cb,'\DeltaF/F (%)');
                title('B2 HM (AUC ↑)','FontWeight','normal'); xlabel(''); ylabel('');
            else
                text(0.5,0.5,'No B2','HorizontalAlignment','center'); axis(ax4,'off');
            end
            axis('square');
            clear o2;
        end
    end
end

%% ===== Helper functions =====
function [alignedData, timeWindow] = extract_aligned_traces(timeVec, signal, tagTimes, pre, post)
    numTags = numel(tagTimes);
    dt = median(diff(timeVec));
    numPts = round((post - pre)/dt) + 1;
    timeWindow = linspace(pre, post, numPts)';  % column
    alignedData = NaN(numTags, numel(timeWindow));
    for k = 1:numTags
        desiredTimes = tagTimes(k) + timeWindow;
        if desiredTimes(1) < timeVec(1) || desiredTimes(end) > timeVec(end)
            continue; % out-of-range trial
        end
        seg = interp1(timeVec, signal, desiredTimes, 'linear', NaN);
        if all(isfinite(seg))
            alignedData(k,:) = seg(:)';
        end
    end
    valid = all(isfinite(alignedData),2);
    alignedData = alignedData(valid,:);
end

function [meanTrace, semTrace] = compute_mean_sem(alignedData)
    meanTrace = mean(alignedData, 1, 'omitnan');
    semTrace  = std(alignedData, 0, 1, 'omitnan') ./ sqrt(size(alignedData,1));
end

function normalizedData = compute_deltaF_over_F_percent(alignedData, timeWindow, preTime, postTime)
    % Returns ΔF/F in PERCENT
    baselineStart = max(preTime, -3);
    baselineEnd   = min(0, postTime);
    idxBase = timeWindow >= baselineStart & timeWindow <= baselineEnd;
    if ~any(idxBase), error('No samples in baseline window.'); end
    F0 = mean(alignedData(:, idxBase), 2);
    nz = F0 < 1e-6;
    if any(nz), F0(nz) = NaN; end
    normalizedData = 100 * (alignedData - F0) ./ F0;  % <-- scale to %
    valid = all(isfinite(normalizedData),2);
    normalizedData = normalizedData(valid,:);
end

function shrink_axes(ax)
    if ~ishandle(ax), return; end
    p = get(ax,'Position'); pad = 0.02;
    p(1) = p(1) + pad;  p(2) = p(2) + pad;
    p(3) = p(3) - 2*pad; p(4) = p(4) - 2*pad;
    set(ax,'Position',p);
end

end
