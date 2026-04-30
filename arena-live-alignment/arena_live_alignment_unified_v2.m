function arena_live_alignment_unified_v2()
% ARENA_LIVE_ALIGNMENT_UNIFIED_V2
% Robust live camera alignment tool for:
%   1) Rectangular arena / open field box
%   2) Circular / elliptical zero-maze arena
%
% Major changes from the original version:
%   - Does NOT directly track only the clicked arena points.
%   - Detects many good visual features around the arena.
%   - Tracks those features and estimates a geometric transform.
%   - Applies the transform to the reference arena geometry.
%   - Zero-maze calibration uses one editable ellipse ROI instead of 10 clicks.
%
% Workflow:
%   - Select webcam
%   - Select arena type
%   - Live preview
%   - Press SPACE to freeze
%   - Calibrate:
%       Rectangle: click 4 inner corners
%       Zero maze: draw one ellipse around the OUTER ring
%   - Live tracking/alignment feedback
%   - Press R to recalibrate
%   - Press ESC to quit
%
% Requirements:
%   - MATLAB Support Package for USB Webcams
%   - Computer Vision Toolbox
%   - Image Processing Toolbox is recommended for drawellipse
%
% Save this file as:
%   arena_live_alignment_unified_v2.m

clc;
close all;

%% ---------------- Parameters ----------------
% Rectangle mode thresholds
rect_angleTol_deg = 1.5;
rect_symTol       = 0.08;
rect_centerTol_px = 20;

% Circle / zero-maze mode thresholds
circ_axisRatioTol = 0.97;   % minor/major should be close to 1 if camera is top-down
circ_centerTol_px = 20;

% Robust tracker parameters
trackingTransformType = 'affine';  % options: 'similarity', 'affine', 'projective'
minValidFeatures     = 12;
refreshInterval      = 60;        % re-seed visual features every N frames
ellipseSampleCount   = 80;        % automatically sampled points on zero-maze ellipse

%% ---------------- Select webcam ----------------
[cam, camName] = selectWebcam();

%% ---------------- Select arena type ----------------
[arenaType, arenaLabel] = selectArenaType();

%% ---------------- Create main figure ----------------
hFig = figure( ...
    'Name', 'Arena Live Alignment Unified V2', ...
    'Color', 'w', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none');

set(hFig, 'KeyPressFcn', @(src,event) keyHandler(src,event));

state.mode = 'preview';              % preview | frozen | tracking
state.freezeRequested = false;
state.stopRequested = false;
state.reselectRequested = false;

state.points = [];
state.refArenaPoints = [];
state.refFeaturePoints = [];
state.tracker = [];

state.camName = camName;
state.arenaType = arenaType;         % 'rectangle' or 'circle'
state.arenaLabel = arenaLabel;
state.frozenFrame = [];
state.lastGoodFrame = [];
state.lastTrackedPoints = [];

state.frameCounter = 0;
state.refreshInterval = refreshInterval;
state.trackingTransformType = trackingTransformType;
state.minValidFeatures = minValidFeatures;
state.ellipseSampleCount = ellipseSampleCount;

setappdata(hFig, 'state', state);

%% ---------------- Main loop ----------------
while ishandle(hFig)
    state = getappdata(hFig, 'state');

    if state.stopRequested
        break;
    end

    switch state.mode

        %% =========================================================
        case 'preview'
            frame = snapshot(cam);

            imshow(frame, 'Border', 'tight');
            hold on;
            drawModeBanner('LIVE PREVIEW MODE', [1 1 0]);
            drawTextLine(60,  'Adjust the camera roughly now.');
            drawTextLine(90,  'Press SPACE to freeze the frame for calibration.');
            drawTextLine(120, 'Press ESC to quit.');
            drawTextLine(150, ['Camera: ' char(state.camName)]);
            drawTextLine(180, ['Arena type: ' char(state.arenaLabel)]);

            if strcmp(state.arenaType, 'rectangle')
                drawTextLine(210, 'Calibration: click 4 INNER arena corners after freeze.');
            else
                drawTextLine(210, 'Calibration: draw one ellipse around the OUTER zero-maze ring after freeze.');
            end
            hold off;
            drawnow limitrate;

            if state.freezeRequested
                state.freezeRequested = false;
                state.frozenFrame = frame;
                state.mode = 'frozen';
                setappdata(hFig, 'state', state);
            end

        %% =========================================================
        case 'frozen'
            frame = state.frozenFrame;

            imshow(frame, 'Border', 'tight');
            hold on;
            drawModeBanner('FROZEN FRAME MODE', [1 1 0]);

            if strcmp(state.arenaType, 'rectangle')
                drawTextLine(60,  'Click 4 INNER arena corners in this order:');
                drawTextLine(90,  'top-left, top-right, bottom-right, bottom-left');
            else
                drawTextLine(60,  'Draw one ellipse around the OUTER zero-maze ring.');
                drawTextLine(90,  'Resize/rotate the ellipse, then double-click inside it or press Enter.');
            end

            drawTextLine(120, 'Press ESC to quit.');
            drawTextLine(150, ['Camera: ' char(state.camName)]);
            drawTextLine(180, ['Arena type: ' char(state.arenaLabel)]);
            hold off;
            drawnow;

            try
                if strcmp(state.arenaType, 'rectangle')
                    points = calibrateRectangleByClicks();
                else
                    points = calibrateZeroMazeByEllipse(gca, state.ellipseSampleCount);
                end
            catch ME
                warning('Calibration cancelled or failed: %s', ME.message);
                state.mode = 'preview';
                setappdata(hFig, 'state', state);
                continue;
            end

            if isempty(points)
                state.mode = 'preview';
                setappdata(hFig, 'state', state);
                continue;
            end

            imshow(frame, 'Border', 'tight');
            hold on;

            if strcmp(state.arenaType, 'rectangle')
                plot([points(:,1); points(1,1)], [points(:,2); points(1,2)], ...
                    'r-o', 'LineWidth', 2, 'MarkerSize', 8);
            else
                [xc, yc, a, b, theta] = fitEllipseToPoints(points);
                drawEllipseOverlay(xc, yc, a, b, theta, 'r', 2);
                plot(xc, yc, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            end

            drawModeBanner('CALIBRATION SELECTED', [0 1 0]);
            drawTextLine(60, 'Detecting robust visual features near the arena...');
            hold off;
            drawnow;

            if ~isempty(state.tracker)
                try
                    release(state.tracker);
                catch
                end
            end

            try
                [tracker, refFeaturePoints] = initializeRobustArenaTracker(frame, points);
            catch ME
                imshow(frame, 'Border', 'tight');
                hold on;
                drawModeBanner('FEATURE DETECTION FAILED', [1 0 0]);
                drawTextLine(60,  'Not enough visual features were detected near the arena.');
                drawTextLine(90,  'Improve lighting/contrast or add visible tape/markers near the arena.');
                drawTextLine(120, 'Press SPACE to try calibration again.');
                drawTextLine(150, ME.message);
                hold off;
                drawnow;

                state.mode = 'preview';
                state.freezeRequested = false;
                setappdata(hFig, 'state', state);
                pause(1.0);
                continue;
            end

            state.points = points;
            state.refArenaPoints = points;
            state.refFeaturePoints = refFeaturePoints;
            state.tracker = tracker;
            state.mode = 'tracking';
            state.frameCounter = 0;
            state.lastTrackedPoints = points;
            state.lastGoodFrame = frame;
            setappdata(hFig, 'state', state);

        %% =========================================================
        case 'tracking'
            frame = snapshot(cam);

            if state.reselectRequested
                state.reselectRequested = false;
                state.freezeRequested = false;
                state.frozenFrame = frame;
                state.mode = 'frozen';
                setappdata(hFig, 'state', state);
                continue;
            end

            if isempty(state.tracker) || isempty(state.refFeaturePoints) || isempty(state.refArenaPoints)
                state.mode = 'preview';
                setappdata(hFig, 'state', state);
                continue;
            end

            try
                [trackedFeaturePoints, validity] = step(state.tracker, frame);
            catch ME
                imshow(frame, 'Border', 'tight');
                hold on;
                drawModeBanner('TRACKER ERROR', [1 0 0]);
                drawTextLine(60,  ME.message);
                drawTextLine(90,  'Press R to recalibrate. Press ESC to quit.');
                hold off;
                drawnow limitrate;
                continue;
            end

            validRef = state.refFeaturePoints(validity, :);
            validCur = trackedFeaturePoints(validity, :);

            imshow(frame, 'Border', 'tight');
            hold on;

            if size(validCur,1) < state.minValidFeatures
                drawModeBanner('TRACKING WEAK', [1 0 0]);
                drawTextLine(60,  sprintf('Only %d reliable features remain.', size(validCur,1)));
                drawTextLine(90,  'Press R to re-freeze and recalibrate.');
                drawTextLine(120, 'Try improving contrast, lighting, or arena-edge visibility.');
                drawTextLine(150, ['Camera: ' char(state.camName)]);
                drawTextLine(180, ['Arena type: ' char(state.arenaLabel)]);
                hold off;
                drawnow limitrate;
                continue;
            end

            try
                tform = estimateGeometricTransform2D( ...
                    validRef, validCur, ...
                    state.trackingTransformType, ...
                    'MaxDistance', 8, ...
                    'Confidence', 95, ...
                    'MaxNumTrials', 200);

                currentArenaPoints = transformPointsForward(tform, state.refArenaPoints);

            catch ME
                drawModeBanner('TRANSFORM FAILED', [1 0 0]);
                drawTextLine(60,  'Feature tracking exists, but transform estimation failed.');
                drawTextLine(90,  'Press R to recalibrate.');
                drawTextLine(120, ME.message);
                hold off;
                drawnow limitrate;
                continue;
            end

            state.lastTrackedPoints = currentArenaPoints;
            state.lastGoodFrame = frame;
            state.frameCounter = state.frameCounter + 1;

            drawModeBanner('LIVE ALIGNMENT MODE', [1 1 0]);

            if strcmp(state.arenaType, 'rectangle')
                corners = reorderCorners(currentArenaPoints);

                metrics = computeRectangleMetrics( ...
                    corners, size(frame), ...
                    rect_angleTol_deg, rect_symTol, rect_centerTol_px);

                feedback = generateRectangleFeedback( ...
                    metrics, rect_angleTol_deg, rect_symTol, rect_centerTol_px);

                plot([corners(:,1); corners(1,1)], [corners(:,2); corners(1,2)], ...
                    'g-o', 'LineWidth', 2, 'MarkerSize', 8);

                imgCenter = [size(frame,2)/2, size(frame,1)/2];
                rectCenter = mean(corners,1);

                plot(imgCenter(1), imgCenter(2), 'bx', 'MarkerSize', 14, 'LineWidth', 2);
                plot(rectCenter(1), rectCenter(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

                drawMetricLine(60,  sprintf('Rotation error: %.2f deg', metrics.rotationDeg));
                drawMetricLine(90,  sprintf('Left-right asymmetry: %.3f', metrics.lrAsym));
                drawMetricLine(120, sprintf('Top-bottom asymmetry: %.3f', metrics.tbAsym));
                drawMetricLine(150, sprintf('Center offset: (%.1f, %.1f) px', ...
                    metrics.centerOffset(1), metrics.centerOffset(2)));
                drawMetricLine(180, sprintf('Valid visual features: %d', size(validCur,1)));
                drawMetricLine(210, ['Camera: ' char(state.camName)]);
                drawMetricLine(240, ['Arena type: ' char(state.arenaLabel)]);

            else
                pts = currentArenaPoints;
                [xc, yc, a, b, theta] = fitEllipseToPoints(pts);

                metrics = computeCircleMetrics( ...
                    xc, yc, a, b, theta, size(frame), ...
                    circ_axisRatioTol, circ_centerTol_px);

                feedback = generateCircleFeedback( ...
                    metrics, circ_axisRatioTol, circ_centerTol_px);

                drawEllipseOverlay(xc, yc, a, b, theta, 'g', 2);

                imgCenter = [size(frame,2)/2, size(frame,1)/2];
                plot(imgCenter(1), imgCenter(2), 'bx', 'MarkerSize', 14, 'LineWidth', 2);
                plot(xc, yc, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

                drawMetricLine(60,  sprintf('Major axis: %.1f px', metrics.majorAxis));
                drawMetricLine(90,  sprintf('Minor axis: %.1f px', metrics.minorAxis));
                drawMetricLine(120, sprintf('Axis ratio minor/major: %.4f', metrics.axisRatio));
                drawMetricLine(150, sprintf('Ellipse angle: %.2f deg', metrics.thetaDeg));
                drawMetricLine(180, sprintf('Center offset: (%.1f, %.1f) px', ...
                    metrics.centerOffset(1), metrics.centerOffset(2)));
                drawMetricLine(210, sprintf('Valid visual features: %d', size(validCur,1)));
                drawMetricLine(240, ['Camera: ' char(state.camName)]);
                drawMetricLine(270, ['Arena type: ' char(state.arenaLabel)]);
            end

            y0 = 320;
            for k = 1:numel(feedback)
                text(20, y0 + 28*(k-1), feedback{k}, ...
                    'Color', 'c', 'FontSize', 13, 'FontWeight', 'bold');
            end

            if metrics.isGood
                rectangle('Position', [20 455 28 28], ...
                    'FaceColor', [0 1 0], 'EdgeColor', 'w');
                text(60, 475, 'Alignment OK', ...
                    'Color', 'g', 'FontSize', 14, 'FontWeight', 'bold');
            else
                rectangle('Position', [20 455 28 28], ...
                    'FaceColor', [1 0 0], 'EdgeColor', 'w');
                text(60, 475, 'Adjust camera', ...
                    'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');
            end

            text(20, 520, 'Press R to re-freeze/recalibrate. Press ESC to quit.', ...
                'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');

            % Periodic feature re-seeding prevents long-term drift.
            if state.frameCounter >= state.refreshInterval
                try
                    release(state.tracker);
                catch
                end

                try
                    [newTracker, newFeaturePoints] = initializeRobustArenaTracker(frame, currentArenaPoints);
                    state.tracker = newTracker;
                    state.refFeaturePoints = newFeaturePoints;
                    state.refArenaPoints = currentArenaPoints;
                    state.frameCounter = 0;
                catch
                    % Keep the old tracker if refresh fails.
                end
            end

            setappdata(hFig, 'state', state);
            hold off;
            drawnow limitrate;
    end
end

%% ---------------- Cleanup ----------------
try
    state = getappdata(hFig, 'state');
    if isfield(state, 'tracker') && ~isempty(state.tracker)
        release(state.tracker);
    end
catch
end

clear cam;

if ishandle(hFig)
    close(hFig);
end

end

%% ========================================================================
function points = calibrateRectangleByClicks()
% Click 4 inner corners. The function reorders them afterward for safety.

nPts = 4;
[x, y, button] = ginput(nPts);

if numel(x) < nPts || numel(y) < nPts
    points = [];
    return;
end

if ~isempty(button) && any(button ~= 1)
    points = [];
    return;
end

points = reorderCorners([x y]);

end

%% ========================================================================
function points = calibrateZeroMazeByEllipse(ax, nSample)
% Preferred method: draw one editable ellipse around the outer ring.
% Fallback method: if drawellipse is not available, click 5-8 points.

if exist('drawellipse', 'file') == 2
    axes(ax);
    hEllipse = drawellipse(ax, ...
        'Color', 'r', ...
        'LineWidth', 2, ...
        'InteractionsAllowed', 'all');

    wait(hEllipse);

    xc = hEllipse.Center(1);
    yc = hEllipse.Center(2);
    a  = hEllipse.SemiAxes(1);
    b  = hEllipse.SemiAxes(2);
    theta = deg2rad(hEllipse.RotationAngle);

    points = sampleEllipsePoints(xc, yc, a, b, theta, nSample);
else
    title('drawellipse not found. Click 6-8 points on the OUTER ring, then press Enter.');
    [x, y] = getpts(ax);
    if numel(x) < 5
        error('At least 5 points are needed for ellipse fallback calibration.');
    end
    rawPts = [x y];
    [xc, yc, a, b, theta] = fitEllipseToPoints(rawPts);
    points = sampleEllipsePoints(xc, yc, a, b, theta, nSample);
end

end

%% ========================================================================
function [tracker, featurePoints] = initializeRobustArenaTracker(frame, arenaPoints)
% Detects strong visual features near the arena and initializes a KLT tracker.
% This is much more stable than directly tracking manually clicked points.

gray = toGray(frame);

xMin = floor(min(arenaPoints(:,1)));
xMax = ceil(max(arenaPoints(:,1)));
yMin = floor(min(arenaPoints(:,2)));
yMax = ceil(max(arenaPoints(:,2)));

margin = 50;

xMin = max(1, xMin - margin);
yMin = max(1, yMin - margin);
xMax = min(size(gray,2), xMax + margin);
yMax = min(size(gray,1), yMax + margin);

roi = [xMin, yMin, max(1, xMax - xMin), max(1, yMax - yMin)];

corners = detectMinEigenFeatures(gray, ...
    'ROI', roi, ...
    'MinQuality', 0.005, ...
    'FilterSize', 5);

if corners.Count < 15
    corners = detectMinEigenFeatures(gray, ...
        'ROI', roi, ...
        'MinQuality', 0.001, ...
        'FilterSize', 3);
end

if corners.Count < 8
    error('Not enough visual features detected near the arena. Add contrast markers or improve lighting.');
end

corners = corners.selectStrongest(min(150, corners.Count));
featurePoints = corners.Location;

tracker = vision.PointTracker( ...
    'MaxBidirectionalError', 8, ...
    'NumPyramidLevels', 4, ...
    'BlockSize', [31 31]);

initialize(tracker, featurePoints, frame);

end

%% ========================================================================
function gray = toGray(frame)

if size(frame,3) == 3
    gray = rgb2gray(frame);
else
    gray = frame;
end

gray = im2uint8(gray);

end

%% ========================================================================
function [cam, selectedCam] = selectWebcam()

cams = webcamlist;

if isempty(cams)
    errordlg('No webcam detected. Please plug in a webcam and try again.', ...
        'Webcam Error');
    error('No webcam detected.');
elseif numel(cams) == 1
    selectedCam = cams{1};
    fprintf('Only one webcam detected. Using: %s\n', selectedCam);
    cam = webcam(selectedCam);
    return;
end

[idx, tf] = listdlg( ...
    'PromptString', 'Select a webcam for calibration:', ...
    'SelectionMode', 'single', ...
    'ListString', cams, ...
    'ListSize', [360 180], ...
    'Name', 'Webcam Selection');

if ~tf
    error('No webcam selected. User cancelled the operation.');
end

selectedCam = cams{idx};

try
    cam = webcam(selectedCam);
    fprintf('Using camera: %s\n', selectedCam);
catch ME
    error('Failed to connect to webcam "%s".\n%s', selectedCam, ME.message);
end

end

%% ========================================================================
function [arenaType, arenaLabel] = selectArenaType()

choices = { ...
    'Rectangular arena (open field / box)', ...
    'Circular / elliptical arena (zero maze) - draw ellipse'};

[idx, tf] = listdlg( ...
    'PromptString', 'Select arena type:', ...
    'SelectionMode', 'single', ...
    'ListString', choices, ...
    'ListSize', [460 140], ...
    'Name', 'Arena Type');

if ~tf
    error('No arena type selected. User cancelled the operation.');
end

if idx == 1
    arenaType = 'rectangle';
    arenaLabel = 'Rectangular arena';
else
    arenaType = 'circle';
    arenaLabel = 'Circular / elliptical zero maze';
end

end

%% ========================================================================
function keyHandler(src, event)

state = getappdata(src, 'state');

switch lower(event.Key)
    case 'space'
        if strcmp(state.mode, 'preview')
            state.freezeRequested = true;
        end
    case 'r'
        if strcmp(state.mode, 'tracking')
            state.reselectRequested = true;
        elseif strcmp(state.mode, 'preview')
            state.freezeRequested = true;
        end
    case 'escape'
        state.stopRequested = true;
end

setappdata(src, 'state', state);

end

%% ========================================================================
function corners = reorderCorners(pts)
% Returns corners in this order:
% top-left, top-right, bottom-right, bottom-left

if size(pts,1) < 4
    corners = pts;
    return;
end

% Use the centroid and angle sorting first.
c = mean(pts, 1);
ang = atan2(pts(:,2) - c(2), pts(:,1) - c(1));
[~, idx] = sort(ang);
ordered = pts(idx, :);

% MATLAB image coordinates have y increasing downward.
% Find top-left as the point with smallest x+y.
s = ordered(:,1) + ordered(:,2);
[~, startIdx] = min(s);
ordered = circshift(ordered, -(startIdx-1), 1);

% Ensure clockwise-like order: TL, TR, BR, BL.
% If second point is not top-right, flip the order after TL.
if ordered(2,1) < ordered(end,1)
    ordered = [ordered(1,:); flipud(ordered(2:end,:))];
end

corners = ordered(1:4, :);

end

%% ========================================================================
function metrics = computeRectangleMetrics(corners, frameSize, angleTol_deg, symTol, centerTol_px)

p1 = corners(1,:);  % top-left
p2 = corners(2,:);  % top-right
p3 = corners(3,:);  % bottom-right
p4 = corners(4,:);  % bottom-left

topVec    = p2 - p1;
bottomVec = p3 - p4;
leftVec   = p4 - p1;
rightVec  = p3 - p2;

topAng = atan2d(topVec(2), topVec(1));
botAng = atan2d(bottomVec(2), bottomVec(1));
rotationDeg = mean([topAng, botAng]);

topLen    = norm(topVec);
bottomLen = norm(bottomVec);
leftLen   = norm(leftVec);
rightLen  = norm(rightVec);

tbAsym = (topLen - bottomLen) / max(topLen, bottomLen);
lrAsym = (leftLen - rightLen) / max(leftLen, rightLen);

rectCenter = mean(corners,1);
imgCenter  = [frameSize(2)/2, frameSize(1)/2];
centerOffset = rectCenter - imgCenter;

metrics.rotationDeg = rotationDeg;
metrics.tbAsym = tbAsym;
metrics.lrAsym = lrAsym;
metrics.centerOffset = centerOffset;

metrics.isGood = ...
    abs(rotationDeg) < angleTol_deg && ...
    abs(tbAsym) < symTol && ...
    abs(lrAsym) < symTol && ...
    abs(centerOffset(1)) < centerTol_px && ...
    abs(centerOffset(2)) < centerTol_px;

end

%% ========================================================================
function feedback = generateRectangleFeedback(m, angleTol_deg, symTol, centerTol_px)

feedback = {};

if m.rotationDeg > angleTol_deg
    feedback{end+1} = 'Rotate camera clockwise slightly';
elseif m.rotationDeg < -angleTol_deg
    feedback{end+1} = 'Rotate camera counterclockwise slightly';
else
    feedback{end+1} = 'Rotation looks good';
end

if m.lrAsym > symTol
    feedback{end+1} = 'Tilt camera slightly to the RIGHT';
elseif m.lrAsym < -symTol
    feedback{end+1} = 'Tilt camera slightly to the LEFT';
else
    feedback{end+1} = 'Left-right tilt looks good';
end

if m.tbAsym > symTol
    feedback{end+1} = 'Tilt camera slightly toward the TOP side';
elseif m.tbAsym < -symTol
    feedback{end+1} = 'Tilt camera slightly toward the BOTTOM side';
else
    feedback{end+1} = 'Top-bottom tilt looks good';
end

if abs(m.centerOffset(1)) > centerTol_px || abs(m.centerOffset(2)) > centerTol_px
    feedback{end+1} = sprintf('Move camera to recenter: dx=%.0f px, dy=%.0f px', ...
        -m.centerOffset(1), -m.centerOffset(2));
else
    feedback{end+1} = 'Centering looks good';
end

end

%% ========================================================================
function metrics = computeCircleMetrics(xc, yc, a, b, theta, frameSize, axisRatioTol, centerTol_px)

majorAxis = max(a, b) * 2;
minorAxis = min(a, b) * 2;
axisRatio = minorAxis / majorAxis;

imgCenter = [frameSize(2)/2, frameSize(1)/2];
centerOffset = [xc yc] - imgCenter;

thetaDeg = rad2deg(theta);

metrics.majorAxis = majorAxis;
metrics.minorAxis = minorAxis;
metrics.axisRatio = axisRatio;
metrics.thetaDeg = thetaDeg;
metrics.centerOffset = centerOffset;

metrics.isGood = ...
    axisRatio >= axisRatioTol && ...
    abs(centerOffset(1)) < centerTol_px && ...
    abs(centerOffset(2)) < centerTol_px;

end

%% ========================================================================
function feedback = generateCircleFeedback(m, axisRatioTol, centerTol_px)

feedback = {};

if m.axisRatio >= axisRatioTol
    feedback{end+1} = 'Circularity looks good';
else
    feedback{end+1} = 'Tilt camera to reduce ellipse distortion';
end

if abs(m.centerOffset(1)) > centerTol_px
    if m.centerOffset(1) > 0
        feedback{end+1} = 'Move camera slightly LEFT to recenter';
    else
        feedback{end+1} = 'Move camera slightly RIGHT to recenter';
    end
else
    feedback{end+1} = 'Left-right centering looks good';
end

if abs(m.centerOffset(2)) > centerTol_px
    if m.centerOffset(2) > 0
        feedback{end+1} = 'Move camera slightly UP to recenter';
    else
        feedback{end+1} = 'Move camera slightly DOWN to recenter';
    end
else
    feedback{end+1} = 'Top-bottom centering looks good';
end

if m.axisRatio < axisRatioTol
    if abs(cosd(m.thetaDeg)) > abs(sind(m.thetaDeg))
        feedback{end+1} = 'Ellipse is stretched mostly left-right; adjust camera tilt on that axis';
    else
        feedback{end+1} = 'Ellipse is stretched mostly top-bottom; adjust camera tilt on that axis';
    end
end

end

%% ========================================================================
function [xc, yc, a, b, theta] = fitEllipseToPoints(pts)
% PCA ellipse approximation.
% This works well for alignment feedback when points are sampled from a
% manually drawn ellipse or transformed from a sampled ellipse.

x = pts(:,1);
y = pts(:,2);

xc = mean(x);
yc = mean(y);

X = [x - xc, y - yc];
C = cov(X);

[V, D] = eig(C);
[dSorted, order] = sort(diag(D), 'descend');
V = V(:, order);
dSorted = dSorted(order);

% For uniformly sampled ellipse boundary points:
% variance along an axis is approximately semiAxis^2 / 2.
a = sqrt(max(dSorted(1), eps)) * sqrt(2);
b = sqrt(max(dSorted(2), eps)) * sqrt(2);

theta = atan2(V(2,1), V(1,1));

if b > a
    tmp = a;
    a = b;
    b = tmp;
    theta = theta + pi/2;
end

end

%% ========================================================================
function pts = sampleEllipsePoints(xc, yc, a, b, theta, nPts)

t = linspace(0, 2*pi, nPts+1);
t(end) = [];

x = a * cos(t);
y = b * sin(t);

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
xy = R * [x; y];

pts = [xy(1,:)' + xc, xy(2,:)' + yc];

end

%% ========================================================================
function drawEllipseOverlay(xc, yc, a, b, theta, colorSpec, lw)

t = linspace(0, 2*pi, 240);
x = a*cos(t);
y = b*sin(t);

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
xy = R * [x; y];

plot(xy(1,:) + xc, xy(2,:) + yc, 'Color', colorSpec, 'LineWidth', lw);

end

%% ========================================================================
function drawModeBanner(txt, colorRGB)

text(20, 30, txt, ...
    'Color', colorRGB, 'FontSize', 16, 'FontWeight', 'bold');

end

%% ========================================================================
function drawTextLine(y, txt)

text(20, y, txt, ...
    'Color', 'c', 'FontSize', 13, 'FontWeight', 'bold');

end

%% ========================================================================
function drawMetricLine(y, txt)

text(20, y, txt, ...
    'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold');

end
