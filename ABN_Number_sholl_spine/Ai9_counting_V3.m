%% Batch Processing Script: AI9+ Cell Counting with Marker-Controlled Watershed
% This script reads a 3-channel TIF, visualizes merged channels,
% segments AI9+ cells, applies a marker-controlled watershed to split touching cells,
% filters by area and shape, and counts the remaining objects.

clear; clc; close all;

%% 1. Input and Parameters
tifImageFile = 'test4.tif';    % Path to your multi-channel TIFF
areaMin      = 20;             % Minimum acceptable area (pixels)
areaMax      = 2000;           % Maximum acceptable area (pixels)
threshold    = 0.7;            % Binarization threshold for AI9 channel
hMin         = 1;              % Suppress shallow minima (controls over-segmentation)
shapeEccMax  = 0.95;           % Maximum eccentricity (filter out non-round objects)

%% 2. Read Multi-Channel TIF
info = imfinfo(tifImageFile);
numPages = numel(info);
if numPages < 3
    error('TIF image does not have at least 3 channels: %s', tifImageFile);
end
img = zeros(info(1).Height, info(1).Width, numPages, 'double');
for k = 1:numPages
    img(:,:,k) = double(imread(tifImageFile, k));
end
DAPI = min(max(img(:,:,1), 0), 5000) / 5000;
Ai9  = min(max(img(:,:,2), 0), 5000) / 5000;

%% 3. Visualize Merged Channels
merged = cat(3, Ai9, zeros(size(Ai9)), DAPI);
figure; imshow(merged);
title('Merged AI9 (red) and DAPI (blue)');

%% 4. Initial Segmentation of AI9+ Cells
mask = imbinarize(Ai9, threshold);
mask = medfilt2(mask, [5,5]);      % Remove noise
mask = imerode(mask, strel('disk',1));
mask = imclearborder(mask);

figure; imshow(mask);
title('Cleaned Binary Mask');

%% 5. Marker-Controlled Watershed for Separation
% Compute distance transform
D = -bwdist(~mask);
% Suppress shallow minima to avoid over-segmentation
D = imhmin(D, hMin);
% Impose minima at background
D(~mask) = -Inf;
% Compute watershed
L = watershed(D);
% Remove ridge lines
mask_w = mask;
mask_w(L == 0) = 0;

figure; imshow(mask_w);
title('Mask after Marker-Controlled Watershed');

%% 6. Connected Components & Filtering
CC = bwconncomp(mask_w, 8);
stats = regionprops(CC, 'Area', 'Centroid', 'Eccentricity');
areas = [stats.Area];
ecc   = [stats.Eccentricity];
% Filter by area and shape
goodIdx = find(areas >= areaMin & areas <= areaMax & ecc <= shapeEccMax);

% Build filtered mask
globalLabels = labelmatrix(CC);
filteredMask = ismember(globalLabels, goodIdx);

% Final stats
tats = regionprops(filteredMask, 'Centroid', 'Area');
centroids = cat(1, tats.Centroid);
filteredAreas = [tats.Area];

%% 7. Display Results
figure; imshow(Ai9); hold on;
plot(centroids(:,1), centroids(:,2), 'r*');
title(sprintf('Filtered AI9+ Cells (N = %d)', numel(centroids)));
hold off;

figure;
histogram(filteredAreas);
xlabel('Area (pixels)'); ylabel('Count');
title('Distribution of Filtered Cell Areas');

%% 8. Output Count
fprintf('Total AI9+ cells (area [%d, %d], ecc <= %.2f): %d\n', ...
    areaMin, areaMax, shapeEccMax, numel(centroids));
