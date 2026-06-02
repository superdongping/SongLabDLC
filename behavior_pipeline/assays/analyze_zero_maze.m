function result = analyze_zero_maze(pairs, outputDir, options)
%ANALYZE_ZERO_MAZE Zero maze open-arm ratio and movement summary.

firstVideo = VideoReader(pairs(1).videoPath);
firstFrame = readFrame(firstVideo);

figure('Name', 'Zero maze ROI calibration', 'NumberTitle', 'off');
imshow(firstFrame);
title('Draw the polygon for the first closed arm ROI, double-click to finish');
closedArmsROI1 = drawpolygon;
wait(closedArmsROI1);
maskROI1 = closedArmsROI1.Position;

title('Draw the polygon for the second closed arm ROI, double-click to finish');
closedArmsROI2 = drawpolygon;
wait(closedArmsROI2);
maskROI2 = closedArmsROI2.Position;
close(gcf);

summaryData = {};
heatmaps = cell(numel(pairs), 1);
trackData = cell(numel(pairs), 1);
maxTimeSpent = 0;

for i = 1:numel(pairs)
    videoInfo = get_video_info(pairs(i).videoPath);
    fps = ceil(videoInfo.frame_rate);
    dlc = load_dlc_csv(pairs(i).csvPath);
    main = get_dlc_bodypart(dlc, options.body_part_main);

    valid = main.likelihood >= options.likelihood_threshold & ...
        limit_tracking_duration(dlc.frames, fps, options.analysis_duration_sec);
    x = main.x(valid);
    y = main.y(valid);
    frames = dlc.frames(valid);

    xEdges = linspace(min(x), max(x), options.heatmap_grid_bins);
    yEdges = linspace(min(y), max(y), options.heatmap_grid_bins);
    timeSpent = histcounts2(x, y, xEdges, yEdges);

    heatmaps{i} = struct('timeSpent', timeSpent, 'xEdges', xEdges, 'yEdges', yEdges);
    trackData{i} = struct('x', x, 'y', y, 'frames', frames, 'fps', fps, 'videoInfo', videoInfo);
    maxTimeSpent = max(maxTimeSpent, max(timeSpent(:)));
end

for i = 1:numel(pairs)
    x = trackData{i}.x;
    y = trackData{i}.y;
    fps = trackData{i}.fps;

    inClosedArm1 = inpolygon(x, y, maskROI1(:, 1), maskROI1(:, 2));
    inClosedArm2 = inpolygon(x, y, maskROI2(:, 1), maskROI2(:, 2));
    inClosedArm = inClosedArm1 | inClosedArm2;

    totalTime = numel(x) / fps;
    timeInClosedArms = sum(inClosedArm) / fps;
    timeInOpenArms = totalTime - timeInClosedArms;
    ratioOpenArms = timeInOpenArms / totalTime;
    totalDistancePixels = sum(hypot(diff(x), diff(y)));

    [~, baseName, ~] = fileparts(pairs(i).videoFile);
    save_zero_maze_qc(pairs(i).videoPath, x, y, maskROI1, maskROI2, outputDir, baseName);
    save_zero_maze_heatmap(heatmaps{i}, maxTimeSpent, options.heatmap_sigma, outputDir, baseName);

    summaryData(end+1, :) = {pairs(i).csvFile, pairs(i).videoFile, totalDistancePixels, ...
        timeInOpenArms, ratioOpenArms, timeInClosedArms}; %#ok<AGROW>
end

summaryTable = cell2table(summaryData, 'VariableNames', ...
    {'CSV_File', 'Video_File', 'Total_Distance_pixels', 'Time_in_Open_Arms_sec', ...
    'Ratio_Open_Arms', 'Time_in_Closed_Arms_sec'});
summaryPath = save_summary_table(summaryTable, outputDir, 'summary_Zero_maze.xlsx');

result = struct();
result.summary_file = summaryPath;
result.closed_arm_roi_1 = maskROI1;
result.closed_arm_roi_2 = maskROI2;
end

function save_zero_maze_qc(videoPath, x, y, maskROI1, maskROI2, outputDir, baseName)
videoObj = VideoReader(videoPath);
firstFrame = readFrame(videoObj);
figure('Name', sprintf('Zero maze trajectory - %s', baseName), 'NumberTitle', 'off');
imshow(firstFrame);
hold on;
scatter(x, y, 8, 'g', 'filled', 'MarkerFaceAlpha', 0.5);
plot([maskROI1(:, 1); maskROI1(1, 1)], [maskROI1(:, 2); maskROI1(1, 2)], 'r-', 'LineWidth', 2);
plot([maskROI2(:, 1); maskROI2(1, 1)], [maskROI2(:, 2); maskROI2(1, 2)], 'b-', 'LineWidth', 2);
title(sprintf('Video: %s - Trajectory', baseName), 'Interpreter', 'none');
hold off;
saveas(gcf, fullfile(outputDir, sprintf('%s_Trajectory_FirstFrame.tif', baseName)));
close(gcf);
end

function save_zero_maze_heatmap(heatmapStruct, maxTimeSpent, sigma, outputDir, baseName)
upperLimit = max(1, 0.8 * maxTimeSpent);
timeSpentNormalized = heatmapStruct.timeSpent / upperLimit;
timeSpentSmoothed = imgaussfilt(timeSpentNormalized, sigma);

figure('Name', sprintf('Zero maze heatmap - %s', baseName), 'NumberTitle', 'off');
imagesc(heatmapStruct.xEdges, heatmapStruct.yEdges, timeSpentSmoothed');
colormap('jet');
colorbar;
clim([0 0.05]);
title(sprintf('Mouse Heatmap - %s', baseName), 'Interpreter', 'none');
axis equal;
set(gca, 'XLim', [min(heatmapStruct.xEdges) max(heatmapStruct.xEdges)], ...
    'YLim', [min(heatmapStruct.yEdges) max(heatmapStruct.yEdges)]);
saveas(gcf, fullfile(outputDir, sprintf('%s_Heatmap_jet_Smoothed.tif', baseName)));
close(gcf);
end
