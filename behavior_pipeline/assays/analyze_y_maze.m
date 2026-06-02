function result = analyze_y_maze(pairs, outputDir, options)
%ANALYZE_Y_MAZE Y-maze arm sequence and spontaneous alternation.

firstVideo = VideoReader(pairs(1).videoPath);
firstFrame = readFrame(firstVideo);

figure('Name', 'Y-maze ROI calibration', 'NumberTitle', 'off');
imshow(firstFrame);

title('Draw polygon for A-arm, then double-click to finish');
roiA = drawpolygon;
wait(roiA);
maskA = roiA.Position;

title('Draw polygon for B-arm, then double-click to finish');
roiB = drawpolygon;
wait(roiB);
maskB = roiB.Position;

title('Draw polygon for C-arm, then double-click to finish');
roiC = drawpolygon;
wait(roiC);
maskC = roiC.Position;
close(gcf);

summaryData = {};

for i = 1:numel(pairs)
    videoInfo = get_video_info(pairs(i).videoPath);
    fps = ceil(videoInfo.frame_rate);
    dlc = load_dlc_csv(pairs(i).csvPath);
    main = get_dlc_bodypart(dlc, options.body_part_main);

    validDuration = limit_tracking_duration(dlc.frames, fps, options.analysis_duration_sec);
    x = main.x(validDuration);
    y = main.y(validDuration);
    likelihood = main.likelihood(validDuration);

    sequence = build_arm_sequence(x, y, likelihood, options.likelihood_threshold, maskA, maskB, maskC);
    alternationPercentage = compute_alternation_percentage(sequence);

    [~, baseName, ~] = fileparts(pairs(i).videoFile);
    save_y_maze_qc(pairs(i).videoPath, x, y, likelihood, options.likelihood_threshold, ...
        maskA, maskB, maskC, outputDir, baseName);

    summaryData(end+1, :) = {pairs(i).csvFile, pairs(i).videoFile, sequence, alternationPercentage}; %#ok<AGROW>
end

summaryTable = cell2table(summaryData, 'VariableNames', ...
    {'CSV_File', 'Video_File', 'Sequence', 'Alternation_Percentage'});
summaryPath = save_summary_table(summaryTable, outputDir, 'Summary_Y_maze.xlsx');

result = struct();
result.summary_file = summaryPath;
result.roi_A = maskA;
result.roi_B = maskB;
result.roi_C = maskC;
end

function sequence = build_arm_sequence(x, y, likelihood, threshold, maskA, maskB, maskC)
sequence = '';
currentState = 'O';

for j = 1:numel(x)
    if likelihood(j) < threshold
        newState = 'O';
    else
        inA = inpolygon(x(j), y(j), maskA(:, 1), maskA(:, 2));
        inB = inpolygon(x(j), y(j), maskB(:, 1), maskB(:, 2));
        inC = inpolygon(x(j), y(j), maskC(:, 1), maskC(:, 2));

        if inA
            newState = 'A';
        elseif inB
            newState = 'B';
        elseif inC
            newState = 'C';
        else
            newState = 'O';
        end
    end

    if newState ~= currentState
        if newState ~= 'O'
            sequence = [sequence newState]; %#ok<AGROW>
        end
        currentState = newState;
    end
end
end

function alternationPercentage = compute_alternation_percentage(sequence)
totalTriplets = length(sequence) - 2;
if totalTriplets <= 0
    alternationPercentage = NaN;
    return;
end

alternationCount = 0;
for j = 1:totalTriplets
    triplet = sequence(j:j+2);
    if numel(unique(triplet)) == 3
        alternationCount = alternationCount + 1;
    end
end
alternationPercentage = (alternationCount / totalTriplets) * 100;
end

function save_y_maze_qc(videoPath, x, y, likelihood, threshold, maskA, maskB, maskC, outputDir, baseName)
videoObj = VideoReader(videoPath);
firstFrame = readFrame(videoObj);
validPoints = likelihood >= threshold;
xTraj = x(validPoints);
yTraj = y(validPoints);

figure('Name', sprintf('Y-maze trajectory - %s', baseName), 'NumberTitle', 'off');
imshow(firstFrame);
hold on;
plot([maskA(:, 1); maskA(1, 1)], [maskA(:, 2); maskA(1, 2)], 'r-', 'LineWidth', 2, 'DisplayName', 'A-arm');
plot([maskB(:, 1); maskB(1, 1)], [maskB(:, 2); maskB(1, 2)], 'b-', 'LineWidth', 2, 'DisplayName', 'B-arm');
plot([maskC(:, 1); maskC(1, 1)], [maskC(:, 2); maskC(1, 2)], 'y-', 'LineWidth', 2, 'DisplayName', 'C-arm');
plot(xTraj, yTraj, 'g-', 'LineWidth', 2, 'DisplayName', 'Trajectory');

if numel(xTraj) > 1
    numArrows = min(20, numel(xTraj) - 1);
    indices = unique(round(linspace(1, numel(xTraj) - 1, numArrows)));
    dx = diff(xTraj);
    dy = diff(yTraj);
    quiver(xTraj(indices), yTraj(indices), dx(indices), dy(indices), 0, ...
        'MaxHeadSize', 2, 'Color', 'g', 'LineWidth', 1, 'AutoScale', 'off');
end

title(sprintf('Video: %s - Trajectory with Y-maze Arms', baseName), 'Interpreter', 'none');
legend('show', 'Location', 'best');
hold off;
saveas(gcf, fullfile(outputDir, sprintf('%s_Trajectory_FirstFrame.tif', baseName)));
close(gcf);
end
