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
scaleLengthPx = hypot(diff(scaleLine.Position(:, 1)), diff(scaleLine.Position(:, 2)));
pixelsPerMeter = scaleLengthPx / options.open_field_scale_m;

centralWidthPixels = options.center_square_m * pixelsPerMeter;
centralArea = [
    openFieldPosition(1) + (openFieldPosition(3) - centralWidthPixels) / 2, ...
    openFieldPosition(2) + (openFieldPosition(4) - centralWidthPixels) / 2, ...
    centralWidthPixels, centralWidthPixels];
close(gcf);

summaryData = {};
heatmapMin = inf;
heatmapMax = -inf;
heatmaps = cell(numel(pairs), 1);
trackData = cell(numel(pairs), 1);

for i = 1:numel(pairs)
    videoInfo = get_video_info(pairs(i).videoPath);
    fps = round(videoInfo.frame_rate);
    dlc = load_dlc_csv(pairs(i).csvPath);
    main = get_dlc_bodypart(dlc, options.body_part_main);

    [x, y, frames] = clean_main_body_trace(main, dlc.frames, fps, options, outerRegionPosition, openFieldPosition);
    trackData{i} = struct('x', x, 'y', y, 'frames', frames, 'fps', fps, 'videoInfo', videoInfo);

    heatmapData = histcounts2(x, y, ...
        linspace(openFieldPosition(1), openFieldPosition(1) + openFieldPosition(3), options.heatmap_grid_size(1)), ...
        linspace(openFieldPosition(2), openFieldPosition(2) + openFieldPosition(4), options.heatmap_grid_size(2)), ...
        'Normalization', 'probability');
    heatmaps{i} = imgaussfilt(heatmapData, options.heatmap_sigma);
    heatmapMin = min(heatmapMin, min(heatmaps{i}(:)));
    heatmapMax = max(heatmapMax, max(heatmaps{i}(:)));
end

adjustedHeatmapMax = heatmapMax * options.heatmap_display_max_fraction;

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
    save_oft_qc(pairs(i).videoPath, x, y, openFieldPosition, centralArea, scaleLine.Position, outputDir, baseName);
    save_oft_heatmap(heatmaps{i}, heatmapMin, adjustedHeatmapMax, outputDir, baseName);

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
result.scale_line_position = scaleLine.Position;
result.pixels_per_meter = pixelsPerMeter;
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

function save_oft_heatmap(heatmapData, heatmapMin, heatmapMax, outputDir, baseName)
figure('Name', sprintf('OFT heatmap - %s', baseName), 'NumberTitle', 'off');
imagesc(flipud(heatmapData), [heatmapMin, heatmapMax]);
axis square;
colormap('jet');
colorbar;
title('Smoothed Heatmap of MainBody Trajectory');
xlabel('X Position');
ylabel('Y Position');
set(gca, 'YDir', 'normal');
saveas(gcf, fullfile(outputDir, sprintf('%s_SmoothedHeatmap_Adjusted.tif', baseName)));
close(gcf);
end
