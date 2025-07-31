% Sholl Analysis GUI Script for Apical Dendrite of Adult-Born DG Neuron
% Dynamic threshold adjustment via slider, manual ROI & soma selection,
% Sholl analysis within ROI, visualization, XLSX & PNG export.

clc; close all; clear all;

%% 1. Select and import the TIFF image
titleDlg = 'Select 2-channel TIFF Image';
[file, path] = uigetfile('*.tif', titleDlg);
if isequal(file,0)
    error('No file selected.');
end
fullFile = fullfile(path, file);
[~, baseName, ~] = fileparts(file);

tifInfo = imfinfo(fullFile);
numPages = numel(tifInfo);
samplesPerPixel = tifInfo(1).SamplesPerPixel;
if numPages >= 2
    ch1 = imread(fullFile, 1);  % DAPI
    ch2 = imread(fullFile, 2);  % TdTomato
elseif samplesPerPixel >= 2
    img = imread(fullFile);
    ch1 = img(:,:,1);
    ch2 = img(:,:,2);
else
    error('The selected TIFF has only one channel or page.');
end

%% 2. ROI and soma selection
figROI = figure('Name','Draw Neuron Boundary');
imshow(ch2, []); title('Draw ROI around one neuron');
roiHandle = imfreehand;
mask = createMask(roiHandle);
close(figROI);

figSoma = figure('Name','Select Soma Center');
imshow(ch2, []); title('Click on neuron soma center');
[x, y] = ginput(1);
cx = x; cy = y;
close(figSoma);

%% 3. Sholl parameters
pixelsPerMicron = 3.2181;
minRadius = 10; maxRadius = 200; stepSize = 10;
radii = minRadius:stepSize:maxRadius;
filterSize = [3,3]; tol = 0.7;
initialThreshold = 0.14;

%% 4. Create GUI
fig = figure('Name','Sholl Analysis GUI','Units','normalized','Position',[0.1 0.1 0.8 0.8]);
axCurve  = subplot(1,2,1,'Parent',fig);
axOverlay = subplot(1,2,2,'Parent',fig);

% Threshold display and slider
lbl = uicontrol('Parent',fig,'Style','text','Units','normalized',...
    'Position',[0.15 0.02 0.1 0.05],'String',sprintf('Threshold: %.2f', initialThreshold));
slider = uicontrol('Parent',fig,'Style','slider','Units','normalized',...
    'Position',[0.3 0.02 0.4 0.05],'Min',0,'Max',1,'Value',initialThreshold);

% Save button
saveBtn = uicontrol('Parent',fig,'Style','pushbutton','Units','normalized',...
    'Position',[0.75 0.02 0.1 0.05],'String','Save','Callback',@saveCallback);

% Store parameters
handles.file = file;
handles.path = path;
handles.baseName = baseName;
handles.ch2 = ch2;
handles.mask = mask;
handles.cx = cx;
handles.cy = cy;
handles.filterSize = filterSize;
handles.radii = radii;
handles.pixelsPerMicron = pixelsPerMicron;
handles.tol = tol;
handles.axCurve = axCurve;
handles.axOverlay = axOverlay;
handles.lbl = lbl;
handles.threshold = initialThreshold;
handles.counts = [];

guidata(fig, handles);

% Initial plotting
guidata(fig, handles);
updatePlots(handles);

% Slider callback
slider.Callback = @sliderCallback;

%% Callback and helper functions
function sliderCallback(src, ~)
    fig = ancestor(src, 'figure');
    h = guidata(fig);
    h.threshold = src.Value;
    h.lbl.String = sprintf('Threshold: %.2f', h.threshold);
    guidata(fig, h);
    updatePlots(h);
end

function updatePlots(h)
    % Apply mask and filter
    img = h.ch2; img(~h.mask) = 0;
    imgMed = medfilt2(img, h.filterSize);

    % Binarize & isolate soma component
    bw = imbinarize(mat2gray(imgMed), h.threshold) & h.mask;
    CC = bwconncomp(bw);
    labels = labelmatrix(CC);
    somaIdx = sub2ind(size(bw), round(h.cy), round(h.cx));
    somaLabel = labels(somaIdx);
    if somaLabel == 0
        warning('No neuron detected at this threshold'); return;
    end
    bwNeuron = (labels == somaLabel);

    % Skeletonize and compute intersections
    skel = bwskel(bwNeuron);
    [X,Y] = meshgrid(1:size(skel,2), 1:size(skel,1));
    distMap = hypot(X - h.cx, Y - h.cy);
    counts = zeros(size(h.radii));
    for ii = 1:numel(h.radii)
        Rpix = h.radii(ii) * h.pixelsPerMicron;
        band = abs(distMap - Rpix) <= h.tol;
        counts(ii) = sum(skel(band));
    end
    h.counts = counts;
    guidata(ancestor(h.axCurve,'figure'), h);

    % Plot curve
    axes(h.axCurve); cla;
    plot(h.radii, counts, '-o','LineWidth',1.5);
    xlabel('Radius (\mum)'); ylabel('Intersections');
    title('Sholl Analysis Curve'); grid on; axis square;

    % Plot overlay
    axes(h.axOverlay); cla;
    imshow(imgMed, []); hold on;
    [ys, xs] = find(skel);
    plot(xs, ys, '.', 'Color',[0 1 0],'MarkerSize',1);
    theta = linspace(0,2*pi,360);
    for jj = 1:numel(h.radii)
        Rp = h.radii(jj)*h.pixelsPerMicron;
        plot(h.cx+Rp*cos(theta), h.cy+Rp*sin(theta),'w:');
    end
    plot(h.cx, h.cy, 'r+','MarkerSize',10,'LineWidth',1.5);
    title('Sholl Overlay'); hold off;
end

function saveCallback(~, ~)
    fig = gcf;
    h = guidata(fig);
    if isempty(h.counts)
        errordlg('Please adjust threshold until neuron is detected before saving.');
        return;
    end
    % Export XLSX
    T = table(h.radii', h.counts', 'VariableNames', {'Radius_um','Intersections'});
    outX = sprintf('Sholl_analysis_%s.xlsx', h.baseName);
    writetable(T, fullfile(h.path, outX));
    disp(['Sholl analysis saved as: ' fullfile(h.path, outX)]);
    % Export overlay PNG
    % Capture overlay axes
    frame = getframe(h.axOverlay);
    outP = sprintf('Sholl_analysis_overlay_of_%s.png', h.baseName);
    imwrite(frame.cdata, fullfile(h.path, outP));
    disp(['Overlay saved as: ' fullfile(h.path, outP)]);
end
