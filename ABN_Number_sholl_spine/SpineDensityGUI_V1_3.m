function SpineDensityGUI_V1_3
clc; clear; close all;
% SpineDensityGUI_V1_3 - Manual dendritic spine counting GUI
% Allows user to draw an ROI, draw the dendrite shaft polyline,
% select spine positions, and export overlay PNG and Excel (spine coords + density summary).

% Constants
pixelsPerMicron = 9.6544;

% 1. File selection
[filename, pathname] = uigetfile('*.tif', 'Select 2-channel TIFF Image');
if isequal(filename,0), return; end
fullpath = fullfile(pathname, filename);
[~, baseName, ~] = fileparts(filename);
% Read channel 2
img = imread(fullpath, 2);

% 2. Create figure with two panels
hFig = figure('Name','Spine Density Manual Counting','NumberTitle','off', ...
    'Units','normalized','Position',[0.2 0.2 0.6 0.6]);
% Left: full image
hAx1 = subplot(1,2,1,'Parent',hFig);
imshow(img, [], 'Parent', hAx1);
title(hAx1,'Raw Image');
hAx1.NextPlot = 'add';
% Right: ROI zoom placeholder
hAx2 = subplot(1,2,2,'Parent',hFig);
axis(hAx2,'off');
title(hAx2,'ROI Zoom & Annotation');
hAx2.NextPlot = 'add';

% 3. UI Controls
uicontrol(hFig,'Style','pushbutton','String','Draw ROI', ...
    'Units','normalized','Position',[0.02 0.02 0.12 0.06],'Callback',@drawROI);
uicontrol(hFig,'Style','pushbutton','String','Draw Shaft', ...
    'Units','normalized','Position',[0.15 0.02 0.12 0.06],'Callback',@drawShaft);
uicontrol(hFig,'Style','pushbutton','String','Select Spines', ...
    'Units','normalized','Position',[0.28 0.02 0.12 0.06],'Callback',@selectSpines);
uicontrol(hFig,'Style','pushbutton','String','Save and Export', ...
    'Units','normalized','Position',[0.41 0.02 0.18 0.06],'Callback',@saveExport);

% 4. Data storage
data.img      = img;
data.pathname = pathname;
data.baseName = baseName;
data.roiPos   = [];
data.roiMask  = [];
data.cropRect = [];
data.shaftPos = [];
data.spineXY  = [];
guidata(hFig,data);

% --- Callback: Draw ROI ---
function drawROI(~,~)
    data = guidata(hFig);
    % Freehand ROI on full image
    roiObj = imfreehand(hAx1);
    roiPos = roiObj.getPosition;
    mask   = createMask(roiObj);
    plot(hAx1, roiPos(:,1), roiPos(:,2), 'y-', 'LineWidth',1);
    % Crop ROI and display in right panel
    minX = min(roiPos(:,1)); maxX = max(roiPos(:,1));
    minY = min(roiPos(:,2)); maxY = max(roiPos(:,2));
    cropRect = [minX, minY, maxX-minX, maxY-minY];
    roiImg = imcrop(data.img, cropRect);
    cla(hAx2); imshow(roiImg, [], 'Parent', hAx2); hold(hAx2,'on');
    axis(hAx2,'image');
    % Store ROI info
    data.roiPos   = roiPos;
    data.roiMask  = mask;
    data.cropRect = cropRect;
    guidata(hFig,data);
end

% --- Callback: Draw Shaft (polyline) ---
function drawShaft(~,~)
    data = guidata(hFig);
    if isempty(data.roiMask)
        errordlg('Draw ROI first','Error'); return;
    end
    % Draw an open polyline for the dendrite shaft
    shaftObj = drawpolyline(hAx2,'Color','r','LineWidth',2);
    % Retrieve polyline coordinates
    shaftPos = shaftObj.Position;  % Nx2 [x,y]
    data.shaftPos = shaftPos;
    % Save back to GUI data
    guidata(hFig,data);
end

% --- Callback: Select Spines ---
function selectSpines(~,~)
    data = guidata(hFig);
    if isempty(data.roiMask)
        errordlg('Draw ROI first','Error'); return;
    end
    [x,y] = getpts(hAx2);
    plot(hAx2, x, y, 'g.', 'MarkerSize',15);
    data.spineXY = [x(:), y(:)];
    guidata(hFig,data);
end

% --- Callback: Save and Export ---
function saveExport(~,~)
    data = guidata(hFig);
    if isempty(data.roiMask) || isempty(data.shaftPos) || isempty(data.spineXY)
        errordlg('Complete ROI, shaft, and spines selection','Error'); return;
    end
    % Compute dendrite length in pixels and microns
    diffs = diff(data.shaftPos);
    lenPx = sum(sqrt(sum(diffs.^2,2)));
    lenUm = lenPx / pixelsPerMicron;
    % Count spines
    numSpines = size(data.spineXY,1);
    density = numSpines / lenUm;
    
    % Save zoom panel PNG
    frame = getframe(hAx2);
    outPNG = fullfile(data.pathname, ['SpineDensity_' data.baseName '.png']);
    imwrite(frame.cdata, outPNG);
    
    % Export spine coordinates (sheet: Coordinates)
    Tcoord = table(data.spineXY(:,1), data.spineXY(:,2), ...
        'VariableNames', {'SpineX_px','SpineY_px'});
    outXLS = fullfile(data.pathname, ['SpineDensity_' data.baseName '.xlsx']);
    writetable(Tcoord, outXLS, 'Sheet', 'Coordinates');
    
    % Export summary (sheet: Summary)
    Tsum = table(numSpines, lenUm, density, ...
        'VariableNames', {'SpineCount','Length_um','Density_spines_per_um'});
    writetable(Tsum, outXLS, 'Sheet', 'Summary');
    
    msgbox(sprintf('Exported:\nPNG: %s\nCoords & Summary: %s', outPNG, outXLS),'Done');
end
end
