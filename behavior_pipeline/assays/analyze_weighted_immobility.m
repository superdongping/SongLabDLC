function result = analyze_weighted_immobility(pairs, outputDir, options, assayLabel, summaryFileName)
%ANALYZE_WEIGHTED_IMMOBILITY Shared FST/TST weighted immobility analysis.

firstVideo = VideoReader(pairs(1).videoPath);
firstFrame = readFrame(firstVideo);

figure('Name', sprintf('%s scaling calibration', assayLabel), 'NumberTitle', 'off');
imshow(firstFrame);
title(sprintf('Draw a line to indicate the reference height (%.2f cm)', options.scale_height_cm));
scaleLine = drawline('Color', 'red', 'LineWidth', 2);
wait(scaleLine);
pixelLength = hypot(diff(scaleLine.Position(:, 1)), diff(scaleLine.Position(:, 2)));
pixelsToCm = options.scale_height_cm / pixelLength;
close(gcf);

summaryData = {};

for i = 1:numel(pairs)
    fprintf('Processing %s file %d/%d: %s\n', assayLabel, i, numel(pairs), pairs(i).csvFile);
    videoInfo = get_video_info(pairs(i).videoPath);
    fps = round(videoInfo.frame_rate);
    dlc = load_dlc_csv(pairs(i).csvPath);

    keep = limit_tracking_duration(dlc.frames, fps, options.analysis_duration_sec);
    frameCount = sum(keep);
    time = (0:frameCount - 1) / fps;

    speeds = compute_bodypart_speeds(dlc, keep, fps, pixelsToCm, options.body_parts);
    speeds = detect_bodypart_immobility(speeds, options.body_parts, options.speed_threshold_cm_s, ...
        options.continuous_duration_sec, fps);
    [immobilityScore, immobileFrames] = compute_weighted_immobility_score(speeds, options.body_parts, ...
        options.weights, options.immobility_threshold_percent);

    [~, baseName, ~] = fileparts(pairs(i).csvFile);
    save_speed_plot(speeds, options.body_parts, time, options.speed_threshold_cm_s, outputDir, baseName);
    save_speed_heatmap(speeds, options.body_parts, time, outputDir, baseName);
    save_immobility_plot(immobilityScore, immobileFrames, time, options.immobility_threshold_percent, ...
        options.analysis_duration_sec, outputDir, baseName);

    analysisEndFrame = min(floor(options.analysis_duration_sec * fps), numel(immobileFrames));
    analysisImmobile = immobileFrames(1:analysisEndFrame);
    totalImmobileTime = sum(analysisImmobile) / fps;

    latencyToFirstImmobility = NaN;
    if any(analysisImmobile)
        firstImmobileFrame = find(analysisImmobile, 1, 'first');
        latencyToFirstImmobility = firstImmobileFrame / fps;
    end

    summaryData(end+1, :) = {pairs(i).csvFile, pairs(i).videoFile, ...
        totalImmobileTime, latencyToFirstImmobility}; %#ok<AGROW>
end

summaryTable = cell2table(summaryData, 'VariableNames', ...
    {'CSV_File', 'Video_File', 'Total_Immobilized_Time_s', 'Latency_to_First_Immobility_s'});
summaryPath = save_summary_table(summaryTable, outputDir, summaryFileName);

result = struct();
result.summary_file = summaryPath;
result.scale_line_position = scaleLine.Position;
result.pixels_to_cm = pixelsToCm;
result.speed_threshold_cm_s = options.speed_threshold_cm_s;
result.immobility_threshold_percent = options.immobility_threshold_percent;
end

function speeds = compute_bodypart_speeds(dlc, keep, fps, pixelsToCm, bodyParts)
speeds = struct();
for bpIdx = 1:numel(bodyParts)
    bodyPart = bodyParts{bpIdx};
    bp = get_dlc_bodypart(dlc, bodyPart);
    xCm = bp.x(keep) * pixelsToCm;
    yCm = bp.y(keep) * pixelsToCm;
    speedTotal = hypot(diff(xCm), diff(yCm)) * fps;

    fieldName = matlab.lang.makeValidName(bodyPart);
    speeds.(fieldName).name = bodyPart;
    speeds.(fieldName).speed_total = speedTotal(:)';
end
end

function speeds = detect_bodypart_immobility(speeds, bodyParts, speedThreshold, continuousDuration, fps)
N = ceil(continuousDuration * fps);
for bpIdx = 1:numel(bodyParts)
    fieldName = matlab.lang.makeValidName(bodyParts{bpIdx});
    speedTotal = speeds.(fieldName).speed_total;
    binary = isfinite(speedTotal) & speedTotal < speedThreshold;
    convResult = conv(double(binary), ones(1, N), 'same');
    speeds.(fieldName).immobilized = convResult >= N;
end
end

function [immobilityScore, immobileFrames] = compute_weighted_immobility_score(speeds, bodyParts, weights, threshold)
numFrames = numel(speeds.(matlab.lang.makeValidName(bodyParts{1})).speed_total);
immobilityScore = zeros(1, numFrames);

for frame = 1:numFrames
    totalScore = 0;
    for bpIdx = 1:numel(bodyParts)
        bodyPart = bodyParts{bpIdx};
        fieldName = matlab.lang.makeValidName(bodyPart);
        if speeds.(fieldName).immobilized(frame)
            totalScore = totalScore + weights.(bodyPart);
        end
    end
    immobilityScore(frame) = totalScore;
end

immobileFrames = immobilityScore >= threshold;
end

function save_speed_plot(speeds, bodyParts, time, speedThreshold, outputDir, baseName)
speedTime = time(2:end);
figure('Name', sprintf('%s - Total Speeds for All Body Parts', baseName), ...
    'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

for bpIdx = 1:numel(bodyParts)
    fieldName = matlab.lang.makeValidName(bodyParts{bpIdx});
    subplot(numel(bodyParts), 1, bpIdx);
    plot(speedTime, speeds.(fieldName).speed_total, 'b', 'LineWidth', 1.5);
    hold on;
    yline(speedThreshold, 'r--', 'LineWidth', 1.5);
    plot(speedTime(speeds.(fieldName).immobilized), ...
        speeds.(fieldName).speed_total(speeds.(fieldName).immobilized), ...
        'ro', 'MarkerSize', 4);
    xlabel('Time (s)');
    ylabel('Speed (cm/s)');
    title(sprintf('%s - Total Speed', bodyParts{bpIdx}), 'Interpreter', 'none');
    grid on;
    hold off;
end

saveas(gcf, fullfile(outputDir, sprintf('%s_Total_Speeds.tiff', baseName)));
close(gcf);
end

function save_speed_heatmap(speeds, bodyParts, time, outputDir, baseName)
speedTime = time(2:end);
speedMatrix = NaN(numel(bodyParts), numel(speedTime));
for bpIdx = 1:numel(bodyParts)
    fieldName = matlab.lang.makeValidName(bodyParts{bpIdx});
    speedMatrix(bpIdx, :) = speeds.(fieldName).speed_total;
end

figure('Name', sprintf('%s - Heatmap of Total Speeds', baseName), ...
    'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);
imagesc(speedTime, 1:numel(bodyParts), speedMatrix);
colorbar;
colormap('hot');
clim([0 100]);
set(gca, 'YTick', 1:numel(bodyParts), 'YTickLabel', bodyParts);
xlabel('Time (s)');
ylabel('Body Parts');
title('Heatmap of Total Speeds (cm/s) for All Body Parts');
grid on;
saveas(gcf, fullfile(outputDir, sprintf('%s_Speed_Heatmap.tiff', baseName)));
close(gcf);
end

function save_immobility_plot(immobilityScore, immobileFrames, time, threshold, analysisDuration, outputDir, baseName)
speedTime = time(2:end);
figure('Name', sprintf('%s - Weighted Immobility Score', baseName), ...
    'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);
plot(speedTime, immobilityScore, 'b', 'LineWidth', 1.5);
hold on;
yline(threshold, 'r--', 'LineWidth', 1.5);
plot(speedTime(immobileFrames), immobilityScore(immobileFrames), 'ro', 'MarkerSize', 4);
xlabel('Time (s)');
ylabel('Immobility Score (%)');
title('Total Weighted Immobility Score Over Time');
grid on;
xlim([0, analysisDuration]);
hold off;
saveas(gcf, fullfile(outputDir, sprintf('%s_Immobility_Score.tiff', baseName)));
close(gcf);
end
