function dlc = load_dlc_csv(csvPath)
%LOAD_DLC_CSV Load a standard DeepLabCut CSV using body-part names.

[~, bodypartsLine, coordsLine] = read_dlc_header_lines(csvPath);
bodyparts = string(strsplit(bodypartsLine, ','));
coords = lower(string(strsplit(coordsLine, ',')));

bodyparts = bodyparts(2:end);
coords = coords(2:end);

data = readmatrix(csvPath, 'NumHeaderLines', 3);
if isempty(data)
    error('No numeric tracking data found in: %s', csvPath);
end

frames = data(:, 1);
uniqueBodyparts = unique(bodyparts, 'stable');
safeNames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(cellstr(uniqueBodyparts)));

dlc = struct();
dlc.csv_path = csvPath;
dlc.frames = frames;
dlc.bodyparts = cellstr(uniqueBodyparts);
dlc.field_names = safeNames;
dlc.bodypart_map = containers.Map(cellstr(uniqueBodyparts), safeNames);

for i = 1:numel(uniqueBodyparts)
    bp = uniqueBodyparts(i);
    idx = find(bodyparts == bp);
    fieldName = safeNames{i};

    dlc.(fieldName) = struct();
    dlc.(fieldName).name = char(bp);
    dlc.(fieldName).x = get_coord(data, idx, coords, 'x');
    dlc.(fieldName).y = get_coord(data, idx, coords, 'y');
    dlc.(fieldName).likelihood = get_coord(data, idx, coords, 'likelihood');
end
end

function [scorerLine, bodypartsLine, coordsLine] = read_dlc_header_lines(csvPath)
fid = fopen(csvPath, 'r');
if fid == -1
    error('Could not open DLC CSV: %s', csvPath);
end
cleanupObj = onCleanup(@() fclose(fid));

scorerLine = fgetl(fid);
bodypartsLine = fgetl(fid);
coordsLine = fgetl(fid);

if ~ischar(scorerLine) || ~ischar(bodypartsLine) || ~ischar(coordsLine)
    error('DLC CSV must contain scorer/bodyparts/coords header rows: %s', csvPath);
end
clear cleanupObj;
end

function values = get_coord(data, idx, coords, coordName)
coordIdx = idx(coords(idx) == coordName);
if isempty(coordIdx)
    values = NaN(size(data, 1), 1);
    return;
end
values = data(:, coordIdx(1) + 1);
end
