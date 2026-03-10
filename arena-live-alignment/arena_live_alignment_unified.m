function arena_live_alignment_unified()
% ARENA_LIVE_ALIGNMENT_UNIFIED
% Unified live camera alignment tool for:
%   1) Rectangular arena (open field box)
%   2) Circular arena (zero maze)
%
% Workflow:
%   - Detect/select webcam
%   - Select arena type
%   - Live preview
%   - SPACE to freeze
%   - Calibrate:
%       Rectangle: click 4 corners
%       Circle:    click 8-12 points on outer edge
%   - Live tracking with dynamic feedback
%   - R to recalibrate
%   - ESC to quit
%
% Requirements:
%   - MATLAB Support Package for USB Webcams
%   - Computer Vision Toolbox

clc;
close all;

%% ---------------- Parameters ----------------
% Rectangle mode thresholds
rect_angleTol_deg = 1.5;
rect_symTol       = 0.08;
rect_centerTol_px = 20;

% Circle mode thresholds
circ_axisRatioTol = 0.97;   % minor/major should be close to 1
circ_centerTol_px = 20;

% Circle mode point count
defaultCirclePoints = 10;   % recommend 8-12 points

%% ---------------- Select webcam ----------------
[cam, camName] = selectWebcam();

%% ---------------- Select arena type ----------------
[arenaType, arenaLabel] = selectArenaType(defaultCirclePoints);

%% ---------------- Create main figure ----------------
hFig = figure( ...
    'Name', 'Arena Live Alignment Unified', ...
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
state.tracker = [];
state.camName = camName;
state.arenaType = arenaType;         % 'rectangle' or 'circle'
state.arenaLabel = arenaLabel;
state.circlePointCount = defaultCirclePoints;
state.frozenFrame = [];
state.lastGoodFrame = [];
state.lastTrackedPoints = [];

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
                drawTextLine(210, 'Calibration step: click 4 inner corners after freeze.');
            else
                drawTextLine(210, sprintf('Calibration step: click %d points on OUTER circular edge after freeze.', ...
                    state.circlePointCount));
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
                drawTextLine(60,  sprintf('Click %d points along the OUTER circular edge.', ...
                    state.circlePointCount));
                drawTextLine(90,  'Spread the clicks around the full circle as evenly as possible.');
            end

            drawTextLine(120, 'Press ESC to quit.');
            drawTextLine(150, ['Camera: ' char(state.camName)]);
            drawTextLine(180, ['Arena type: ' char(state.arenaLabel)]);
            hold off;
            drawnow;

            try
                if strcmp(state.arenaType, 'rectangle')
                    nPts = 4;
                else
                    nPts = state.circlePointCount;
                end

                [x, y, button] = ginput(nPts);
            catch
                break;
            end

            if numel(x) < nPts || numel(y) < nPts
                state.mode = 'preview';
                setappdata(hFig, 'state', state);
                continue;
            end

            if ~isempty(button) && any(button ~= 1)
                state.mode = 'preview';
                setappdata(hFig, 'state', state);
                continue;
            end

            points = [x y];

            if strcmp(state.arenaType, 'rectangle')
                points = reorderCorners(points);
            end

            imshow(frame, 'Border', 'tight');
            hold on;

            if strcmp(state.arenaType, 'rectangle')
                plot([points(:,1); points(1,1)], [points(:,2); points(1,2)], ...
                    'r-o', 'LineWidth', 2, 'MarkerSize', 8);
            else
                plot(points(:,1), points(:,2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
                [xc, yc, a, b, theta] = fitEllipseToPoints(points);
                drawEllipseOverlay(xc, yc, a, b, theta, 'r', 2);
            end

            drawModeBanner('CALIBRATION POINTS SELECTED', [0 1 0]);
            drawTextLine(60, 'Initializing tracker...');
            hold off;
            drawnow;

            if ~isempty(state.tracker)
                try
                    release(state.tracker);
                catch
                end
            end

            tracker = vision.PointTracker('MaxBidirectionalError', 2);
            initialize(tracker, points, frame);

            state.points = points;
            state.tracker = tracker;
            state.mode = 'tracking';
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

            [trackedPoints, validity] = step(state.tracker, frame);

            imshow(frame, 'Border', 'tight');
            hold on;

            if sum(validity) < size(state.points,1)
                drawModeBanner('TRACKING LOST', [1 0 0]);
                drawTextLine(60,  'Tracking failed.');
                drawTextLine(90,  'Press R to re-freeze and recalibrate.');
                drawTextLine(120, 'Press ESC to quit.');
                drawTextLine(150, ['Camera: ' char(state.camName)]);
                drawTextLine(180, ['Arena type: ' char(state.arenaLabel)]);
                hold off;
                drawnow limitrate;
                continue;
            end

            state.lastTrackedPoints = trackedPoints;
            state.lastGoodFrame = frame;
            setappdata(hFig, 'state', state);

            drawModeBanner('LIVE ALIGNMENT MODE', [1 1 0]);

            if strcmp(state.arenaType, 'rectangle')
                corners = reorderCorners(trackedPoints);

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
                drawMetricLine(180, ['Camera: ' char(state.camName)]);
                drawMetricLine(210, ['Arena type: ' char(state.arenaLabel)]);

            else
                pts = trackedPoints;
                [xc, yc, a, b, theta] = fitEllipseToPoints(pts);

                metrics = computeCircleMetrics( ...
                    xc, yc, a, b, theta, size(frame), ...
                    circ_axisRatioTol, circ_centerTol_px);

                feedback = generateCircleFeedback( ...
                    metrics, circ_axisRatioTol, circ_centerTol_px);

                plot(pts(:,1), pts(:,2), 'go', 'MarkerSize', 7, 'LineWidth', 1.5);
                drawEllipseOverlay(xc, yc, a, b, theta, 'g', 2);

                imgCenter = [size(frame,2)/2, size(frame,1)/2];
                plot(imgCenter(1), imgCenter(2), 'bx', 'MarkerSize', 14, 'LineWidth', 2);
                plot(xc, yc, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

                drawMetricLine(60,  sprintf('Major axis: %.1f px', metrics.majorAxis));
                drawMetricLine(90,  sprintf('Minor axis: %.1f px', metrics.minorAxis));
                drawMetricLine(120, sprintf('Axis ratio (minor/major): %.4f', metrics.axisRatio));
                drawMetricLine(150, sprintf('Ellipse angle: %.2f deg', metrics.thetaDeg));
                drawMetricLine(180, sprintf('Center offset: (%.1f, %.1f) px', ...
                    metrics.centerOffset(1), metrics.centerOffset(2)));
                drawMetricLine(210, ['Camera: ' char(state.camName)]);
                drawMetricLine(240, ['Arena type: ' char(state.arenaLabel)]);
            end

            y0 = 290;
            for k = 1:numel(feedback)
                text(20, y0 + 28*(k-1), feedback{k}, ...
                    'Color', 'c', 'FontSize', 13, 'FontWeight', 'bold');
            end

            if metrics.isGood
                rectangle('Position', [20 430 28 28], ...
                    'FaceColor', [0 1 0], 'EdgeColor', 'w');
                text(60, 450, 'Alignment OK', ...
                    'Color', 'g', 'FontSize', 14, 'FontWeight', 'bold');
            else
                rectangle('Position', [20 430 28 28], ...
                    'FaceColor', [1 0 0], 'EdgeColor', 'w');
                text(60, 450, 'Adjust camera', ...
                    'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');
            end

            text(20, 490, 'Press R to re-freeze / recalibrate. Press ESC to quit.', ...
                'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold');

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
function [arenaType, arenaLabel] = selectArenaType(defaultCirclePoints)
choices = { ...
    'Rectangular arena (open field / box)', ...
    sprintf('Circular arena (zero maze) - click %d edge points', defaultCirclePoints)};

[idx, tf] = listdlg( ...
    'PromptString', 'Select arena type:', ...
    'SelectionMode', 'single', ...
    'ListString', choices, ...
    'ListSize', [420 140], ...
    'Name', 'Arena Type');

if ~tf
    error('No arena type selected. User cancelled the operation.');
end

if idx == 1
    arenaType = 'rectangle';
    arenaLabel = 'Rectangular arena';
else
    arenaType = 'circle';
    arenaLabel = 'Circular arena';
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
% top-left, top-right, bottom-right, bottom-left
s = pts(:,1) + pts(:,2);
d = pts(:,1) - pts(:,2);

topLeft     = pts(find(s == min(s), 1), :);
bottomRight = pts(find(s == max(s), 1), :);
topRight    = pts(find(d == max(d), 1), :);
bottomLeft  = pts(find(d == min(d), 1), :);

corners = [topLeft; topRight; bottomRight; bottomLeft];
end

%% ========================================================================
function metrics = computeRectangleMetrics(corners, frameSize, angleTol_deg, symTol, centerTol_px)

p1 = corners(1,:);
p2 = corners(2,:);
p3 = corners(3,:);
p4 = corners(4,:);

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
% Fit ellipse using covariance/PCA approximation.
% Good for alignment guidance from user-clicked edge points.

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

% Semi-axis estimates from point spread
a = sqrt(max(dSorted(1), eps)) * sqrt(2);
b = sqrt(max(dSorted(2), eps)) * sqrt(2);

theta = atan2(V(2,1), V(1,1));

% Ensure a >= b
if b > a
    tmp = a; a = b; b = tmp;
    theta = theta + pi/2;
end
end

%% ========================================================================
function drawEllipseOverlay(xc, yc, a, b, theta, colorSpec, lw)
t = linspace(0, 2*pi, 200);
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