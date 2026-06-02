function result = analyze_npr(pairs, outputDir, options)
%ANALYZE_NPR Novel object recognition / preference analysis.

firstVideo = VideoReader(pairs(1).videoPath);
firstFrame = readFrame(firstVideo);

figure('Name', 'NPR ROI calibration', 'NumberTitle', 'off');
imshow(firstFrame);
title('Draw a rectangle to define the open field area');
openFieldRect = drawrectangle('Color', 'g');
wait(openFieldRect);
openFieldPosition = openFieldRect.Position;

title(sprintf('Draw a line representing %.2f m', options.open_field_scale_m));
rulerLine = drawline('Color', 'y');
wait(rulerLine);
rulerPixelLength = hypot(diff(rulerLine.Position(:, 1)), diff(rulerLine.Position(:, 2)));
pixelsPerMeter = rulerPixelLength / options.open_field_scale_m;
roiRadiusPixels = options.object_roi_radius_m * pixelsPerMeter;

title('Click on the center of the first object ROI');
[center1_x, center1_y] = ginput(1);
center1 = [center1_x, center1_y];

title('Click on the center of the second object ROI');
[center2_x, center2_y] = ginput(1);
center2 = [center2_x, center2_y];
close(gcf);

summaryData = {};
eventLogData = {};

for i = 1:numel(pairs)
    videoInfo = get_video_info(pairs(i).videoPath);
    fps = videoInfo.frame_rate;
    dlc = load_dlc_csv(pairs(i).csvPath);

    nose = get_dlc_bodypart(dlc, options.body_part_nose);
    head = get_dlc_bodypart(dlc, options.body_part_head);
    body = get_dlc_bodypart(dlc, options.body_part_main);

    valid = nose.likelihood >= options.likelihood_threshold & ...
        head.likelihood >= options.likelihood_threshold & ...
        body.likelihood >= options.likelihood_threshold & ...
        limit_tracking_duration(dlc.frames, fps, options.analysis_duration_sec);

    noseX = nose.x(valid);
    noseY = nose.y(valid);
    bodyX = body.x(valid);
    bodyY = body.y(valid);
    framesValid = dlc.frames(valid);

    inOpen = noseX >= openFieldPosition(1) & noseX <= (openFieldPosition(1) + openFieldPosition(3)) & ...
        noseY >= openFieldPosition(2) & noseY <= (openFieldPosition(2) + openFieldPosition(4));
    noseX = noseX(inOpen);
    noseY = noseY(inOpen);
    bodyX = bodyX(inOpen);
    bodyY = bodyY(inOpen);
    framesValid = framesValid(inOpen);

    inROI1 = hypot(noseX - center1(1), noseY - center1(2)) <= roiRadiusPixels;
    inROI2 = hypot(noseX - center2(1), noseY - center2(2)) <= roiRadiusPixels;

    [exploreROI1, eventsROI1, totalTimeROI1, timestampsROI1] = detect_exploration_bouts( ...
        inROI1, framesValid, center1, roiRadiusPixels, bodyX, bodyY, noseX, noseY, fps, options);
    [exploreROI2, eventsROI2, totalTimeROI2, timestampsROI2] = detect_exploration_bouts( ...
        inROI2, framesValid, center2, roiRadiusPixels, bodyX, bodyY, noseX, noseY, fps, options);

    [~, baseName, ~] = fileparts(pairs(i).videoFile);
    save_npr_qc(pairs(i).videoPath, noseX, noseY, exploreROI1, exploreROI2, ...
        center1, center2, roiRadiusPixels, openFieldPosition, outputDir, baseName);

    summaryData(end+1, :) = {pairs(i).csvFile, pairs(i).videoFile, ...
        totalTimeROI1, eventsROI1, totalTimeROI2, eventsROI2}; %#ok<AGROW>
    eventLogData(end+1, :) = {pairs(i).csvFile, timestampsROI1, timestampsROI2}; %#ok<AGROW>
end

summaryTable = cell2table(summaryData, 'VariableNames', ...
    {'CSV_File', 'Video_File', 'Time_Object_1_sec', 'Events_Object_1', 'Time_Object_2_sec', 'Events_Object_2'});
summaryTable = add_npr_ratios_and_notes(summaryTable);
summaryPath = save_summary_table(summaryTable, outputDir, 'summary_NPR.xlsx', 'Sheet', 'Summary');

eventLogTable = build_event_log_table(eventLogData);
writetable(eventLogTable, summaryPath, 'Sheet', 'Event_Logs');

result = struct();
result.summary_file = summaryPath;
result.open_field_position = openFieldPosition;
result.ruler_line_position = rulerLine.Position;
result.pixels_per_meter = pixelsPerMeter;
result.object_roi_radius_pixels = roiRadiusPixels;
result.object_roi_centers = [center1; center2];
end

function [exploreFlags, numEvents, totalTime, timestamps] = detect_exploration_bouts( ...
    inROI, framesValid, center, radius, bodyX, bodyY, noseX, noseY, fps, options)

N = numel(inROI);
exploreFlags = false(N, 1);
timestamps = [];
if N == 0
    numEvents = 0;
    totalTime = 0;
    return;
end

windowSize = max(1, round(options.window_duration_sec * fps));
minFramesInROI = ceil(options.min_overlap * windowSize);
roiCount = movsum(double(inROI), [windowSize - 1, 0]);
candidate = roiCount >= minFramesInROI;

dist = hypot(noseX - center(1), noseY - center(2));
dDist = [0; diff(dist)];
approach = movmean(dDist, [windowSize - 1, 0]) < 0;
candidate = candidate & approach;

stableFrames = max(1, round(options.sitting_stability_sec * fps));
bodyDist = hypot(bodyX - center(1), bodyY - center(2));
overlap = bodyDist <= (radius / 2);
stable = movstd(bodyX, [stableFrames - 1, 0], 'omitnan') < options.movement_threshold_px & ...
    movstd(bodyY, [stableFrames - 1, 0], 'omitnan') < options.movement_threshold_px;

exploreRaw = candidate & ~(overlap & stable);
exploreMerged = merge_logical_bouts(exploreRaw, max(0, round(options.merge_gap_sec * fps)));
exploreFlags = remove_short_bouts(exploreMerged, max(1, round(options.min_bout_sec * fps)));

[starts, ends] = logical_bouts(exploreFlags);
numEvents = numel(starts);
if numEvents == 0
    totalTime = 0;
    return;
end

durationsSec = (ends - starts + 1) / fps;
totalTime = sum(durationsSec);
timestamps = framesValid(starts) / fps;
end

function out = merge_logical_bouts(x, maxGap)
out = x(:);
if maxGap <= 0
    return;
end
[s, e] = logical_bouts(out);
for k = 1:numel(s) - 1
    gap = s(k + 1) - e(k) - 1;
    if gap > 0 && gap <= maxGap
        out(e(k) + 1:s(k + 1) - 1) = true;
    end
end
end

function out = remove_short_bouts(x, minLen)
out = x(:);
[s, e] = logical_bouts(out);
for k = 1:numel(s)
    if (e(k) - s(k) + 1) < minLen
        out(s(k):e(k)) = false;
    end
end
end

function [starts, ends] = logical_bouts(x)
x = x(:);
d = diff([false; x; false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
end

function save_npr_qc(videoPath, noseX, noseY, exploreROI1, exploreROI2, center1, center2, radius, openFieldPosition, outputDir, baseName)
videoObj = VideoReader(videoPath);
firstFrame = readFrame(videoObj);
figure('Name', sprintf('NPR exploration trajectory - %s', baseName), 'NumberTitle', 'off');
imshow(firstFrame);
hold on;
plot(noseX, noseY, 'm', 'LineWidth', 1);
plot(noseX(exploreROI1), noseY(exploreROI1), 'b.', 'MarkerSize', 10);
plot(noseX(exploreROI2), noseY(exploreROI2), 'c.', 'MarkerSize', 10);
viscircles(center1, radius, 'EdgeColor', 'r');
viscircles(center2, radius, 'EdgeColor', 'r');
rectangle('Position', openFieldPosition, 'EdgeColor', 'g', 'LineWidth', 2);
title(sprintf('Video: %s', baseName), 'Interpreter', 'none');
hold off;
saveas(gcf, fullfile(outputDir, sprintf('%s_FirstFrame_with_ROIs_and_Exploration.tif', baseName)));
close(gcf);
end

function summaryTable = add_npr_ratios_and_notes(summaryTable)
t1 = summaryTable.Time_Object_1_sec;
t2 = summaryTable.Time_Object_2_sec;
den = t1 + t2;

summaryTable.DR_Object1_vs_Object2 = NaN(height(summaryTable), 1);
summaryTable.DR_Object2_vs_Object1 = NaN(height(summaryTable), 1);
validDen = den > 0;
summaryTable.DR_Object1_vs_Object2(validDen) = (t1(validDen) - t2(validDen)) ./ den(validDen);
summaryTable.DR_Object2_vs_Object1(validDen) = (t2(validDen) - t1(validDen)) ./ den(validDen);

notes = strings(height(summaryTable), 1);
zeroMask = (t1 == 0) | (t2 == 0);
notes(zeroMask) = notes(zeroMask) + "WARNING: zero exploring time detected. ";
lowTotalMask = den < 5;
notes(lowTotalMask) = notes(lowTotalMask) + "WARNING: total exploration time < 5 s. ";
denZeroMask = den == 0;
notes(denZeroMask) = notes(denZeroMask) + "WARNING: discrimination ratio undefined. ";
summaryTable.Note = strtrim(notes);
end

function eventLogTable = build_event_log_table(eventLogData)
csvFile = eventLogData(:, 1);
roi1 = cell(size(csvFile));
roi2 = cell(size(csvFile));
for i = 1:numel(csvFile)
    roi1{i} = format_timestamps(eventLogData{i, 2});
    roi2{i} = format_timestamps(eventLogData{i, 3});
end
eventLogTable = table(csvFile, roi1, roi2, 'VariableNames', ...
    {'CSV_File', 'Exploration_ROI1_Timestamps', 'Exploration_ROI2_Timestamps'});
end

function out = format_timestamps(values)
if isempty(values)
    out = '';
    return;
end
parts = arrayfun(@(v) sprintf('%.2f', v), values, 'UniformOutput', false);
out = strjoin(parts, ', ');
end
