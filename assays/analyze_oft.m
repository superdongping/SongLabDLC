function result = analyze_oft(pairs, outputDir, options)
%ANALYZE_OFT Open Field Test analysis using standardized DLC loading.

firstVideo = VideoReader(pairs(1).videoPath);
firstFrame = readFrame(firstVideo);

figure('Name', 'OFT ROI calibration', 'NumberTitle', 'off');
imshow(firstFrame);
title('Draw a rectangle to define the overall region for analysis');
outerRegionRect = drawrectangle('Color', 'c');
wait(outerRegionRect);
outerRegionPosition = outerRegionRect.Position;

title('Draw a rectangle to define the open field area');
openFieldRect = drawrectangle('Color', 'g');
wait(openFieldRect);
openFieldPosition = openFieldRect.Position;

title(sprintf('Draw a line to define %.2f m scale on one side of the open field', options.open_field_scale_m));
scaleLine = drawline('Color', 'b');
wait(scaleLine);
scaleLinePosition = scaleLine.Position;
scaleLengthPx = hypot(diff(scaleLinePosition(:, 1)), diff(scaleLinePosition(:, 2)));
pixelsPerMeter = scaleLengthPx / options.open_field_scale_m;

centralWidthPixels = options.center_square_m * pixelsPerMeter;
centralArea = [
    openFieldPosition(1) + (openFieldPosition(3) - centralWidthPixels) / 2, ...
    openFieldPosition(2) + (openFieldPosition(4) - centralWidthPixels) / 2, ...
    centralWidthPixels, centralWidthPixels];
close(gcf);

summaryData = {};
heatmapMax = -inf;
heatmaps = cell(numel(pairs), 1);
trackData = cell(numel(pairs), 1);
allHeatmapValues = [];
pixelsToCm = 100 / pixelsPerMeter;
openFieldWidthCm = openFieldPosition(3) * pixelsToCm;
openFieldHeightCm = openFieldPosition(4) * pixelsToCm;
xEdgesCm = linspace(0, openFieldWidthCm, options.heatmap_grid_size(1));
yEdgesCm = linspace(0, openFieldHeightCm, options.heatmap_grid_size(2));
xCentersCm = edge_centers(xEdgesCm);
yCentersCm = edge_centers(yEdgesCm);

for i = 1:numel(pairs)
    videoInfo = get_video_info(pairs(i).videoPath);
    fps = round(videoInfo.frame_rate);
    dlc = load_dlc_csv(pairs(i).csvPath);
    main = get_dlc_bodypart(dlc, options.body_part_main);

    [x, y, frames] = clean_main_body_trace(main, dlc.frames, fps, options, outerRegionPosition, openFieldPosition);
    trackData{i} = struct('x', x, 'y', y, 'frames', frames, 'fps', fps, 'videoInfo', videoInfo);

    xCm = (x - openFieldPosition(1)) * pixelsToCm;
    % Video pixel Y increases downward; convert to cm from the arena bottom
    % so the heatmap has the same visual orientation as the trajectory image.
    yCm = (openFieldPosition(2) + openFieldPosition(4) - y) * pixelsToCm;
    occupancyFrames = histcounts2(xCm, yCm, xEdgesCm, yEdgesCm);
    occupancySeconds = occupancyFrames / fps;
    heatmaps{i} = imgaussfilt(occupancySeconds, options.heatmap_sigma);
    heatmapMax = max(heatmapMax, max(heatmaps{i}(:)));
    allHeatmapValues = [allHeatmapValues; heatmaps{i}(:)]; %#ok<AGROW>
end

if isempty(options.heatmap_color_limit_sec)
    heatmapColorMax = robust_heatmap_limit(allHeatmapValues, options.heatmap_color_limit_percentile, heatmapMax);
else
    heatmapColorMax = options.heatmap_color_limit_sec;
end

for i = 1:numel(pairs)
    x = trackData{i}.x;
    y = trackData{i}.y;
    fps = trackData{i}.fps;

    trajectoryDistanceM = sum(hypot(diff(x), diff(y))) / pixelsPerMeter;
    inCentralArea = x >= centralArea(1) & x <= (centralArea(1) + centralArea(3)) & ...
        y >= centralArea(2) & y <= (centralArea(2) + centralArea(4));

    timeInCentralArea = sum(inCentralArea) / fps;
    totalTime = numel(x) / fps;
    centralTimeRatio = timeInCentralArea / totalTime;

    centralX = x(inCentralArea);
    centralY = y(inCentralArea);
    distanceInCentralM = sum(hypot(diff(centralX), diff(centralY))) / pixelsPerMeter;
    speedTotal = trajectoryDistanceM / totalTime;
    speedCentral = distanceInCentralM / timeInCentralArea;

    [~, baseName, ~] = fileparts(pairs(i).videoFile);
    save_oft_qc(pairs(i).videoPath, x, y, openFieldPosition, centralArea, scaleLinePosition, outputDir, baseName);
    save_oft_heatmap(heatmaps{i}, xCentersCm, yCentersCm, heatmapColorMax, outputDir, baseName);

    summaryData(end+1, :) = {pairs(i).csvFile, pairs(i).videoFile, trajectoryDistanceM, ...
        timeInCentralArea, centralTimeRatio, distanceInCentralM, speedTotal, speedCentral}; %#ok<AGROW>
end

summaryTable = cell2table(summaryData, 'VariableNames', {'CSV_File', 'Video_File', ...
    'Total_Trajectory_Distance_m', 'Time_in_Central_s', 'Central_Time_Ratio', ...
    'Distance_in_Central_m', 'Speed_Total_mps', 'Speed_Central_mps'});
summaryPath = save_summary_table(summaryTable, outputDir, 'summary_Open_field.xlsx');

result = struct();
result.summary_file = summaryPath;
result.open_field_position = openFieldPosition;
result.outer_region_position = outerRegionPosition;
result.central_area = centralArea;
result.scale_line_position = scaleLinePosition;
result.pixels_per_meter = pixelsPerMeter;
result.pixels_to_cm = pixelsToCm;
result.open_field_width_cm = openFieldWidthCm;
result.open_field_height_cm = openFieldHeightCm;
result.heatmap_color_limit_sec = heatmapColorMax;
result.heatmap_raw_max_sec = heatmapMax;
result.heatmap_color_limit_percentile = options.heatmap_color_limit_percentile;
end

function [x, y, frames] = clean_main_body_trace(main, frames, fps, options, outerRegionPosition, openFieldPosition)
valid = main.likelihood >= options.likelihood_threshold;
valid = valid & limit_tracking_duration(frames, fps, options.analysis_duration_sec);

x = main.x(valid);
y = main.y(valid);
frames = frames(valid);

inOuter = x >= outerRegionPosition(1) & x <= (outerRegionPosition(1) + outerRegionPosition(3)) & ...
    y >= outerRegionPosition(2) & y <= (outerRegionPosition(2) + outerRegionPosition(4));
x = x(inOuter);
y = y(inOuter);
frames = frames(inOuter);

inOpen = x >= openFieldPosition(1) & x <= (openFieldPosition(1) + openFieldPosition(3)) & ...
    y >= openFieldPosition(2) & y <= (openFieldPosition(2) + openFieldPosition(4));
x = x(inOpen);
y = y(inOpen);
frames = frames(inOpen);
end

function save_oft_qc(videoPath, x, y, openFieldPosition, centralArea, scaleLinePosition, outputDir, baseName)
videoObj = VideoReader(videoPath);
firstFrame = readFrame(videoObj);
figure('Name', sprintf('OFT trajectory - %s', baseName), 'NumberTitle', 'off');
imshow(firstFrame);
hold on;
plot(x, y, 'm', 'LineWidth', 1);
rectangle('Position', openFieldPosition, 'EdgeColor', 'g', 'LineWidth', 2);
rectangle('Position', centralArea, 'EdgeColor', 'r', 'LineWidth', 2);
line(scaleLinePosition(:, 1), scaleLinePosition(:, 2), 'Color', 'b', 'LineWidth', 2);
title(sprintf('Video: %s', baseName), 'Interpreter', 'none');
hold off;
saveas(gcf, fullfile(outputDir, sprintf('%s_FirstFrame_with_Trajectory.tif', baseName)));
close(gcf);
end

function centers = edge_centers(edges)
centers = edges(1:end-1) + diff(edges) / 2;
end

function limitValue = robust_heatmap_limit(values, percentileValue, fallbackMax)
values = values(isfinite(values) & values > 0);
if isempty(values) || isempty(percentileValue) || isnan(percentileValue)
    limitValue = fallbackMax;
    return;
end

values = sort(values(:));
idx = max(1, min(numel(values), ceil((percentileValue / 100) * numel(values))));
limitValue = values(idx);
if limitValue <= 0
    limitValue = fallbackMax;
end
end

function save_oft_heatmap(heatmapData, xCentersCm, yCentersCm, heatmapMax, outputDir, baseName)
figure('Name', sprintf('OFT heatmap - %s', baseName), 'NumberTitle', 'off');
imagesc(xCentersCm, yCentersCm, heatmapData', [0, heatmapMax]);
axis square;
colormap('jet');
c = colorbar;
c.Label.String = sprintf('Time spent (s/bin), capped at %.2g s', heatmapMax);
title('Smoothed Occupancy Heatmap of MainBody Trajectory');
xlabel('X Position (cm)');
ylabel('Y Position (cm from bottom)');
set(gca, 'YDir', 'normal');
saveas(gcf, fullfile(outputDir, sprintf('%s_SmoothedOccupancyHeatmap_sec.tif', baseName)));
close(gcf);
end
