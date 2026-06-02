function pairs = match_csv_video_files(dataFolder)
%MATCH_CSV_VIDEO_FILES Match DLC CSV files to MP4 videos in a folder.

csvFiles = dir(fullfile(dataFolder, '*.csv'));
videoFiles = dir(fullfile(dataFolder, '*.mp4'));

csvFiles = csvFiles(~startsWith({csvFiles.name}, 'summary_', 'IgnoreCase', true));

if isempty(csvFiles)
    error('No CSV files found in: %s', dataFolder);
end
if isempty(videoFiles)
    error('No MP4 files found in: %s', dataFolder);
end

pairs = struct('csvFile', {}, 'videoFile', {}, 'csvPath', {}, ...
    'videoPath', {}, 'baseName', {}, 'matchMethod', {});

usedVideos = false(numel(videoFiles), 1);

for i = 1:numel(csvFiles)
    csvBase = strip_extension(csvFiles(i).name);
    csvNorm = normalize_name(csvBase);
    matchIdx = [];
    matchMethod = '';

    for j = 1:numel(videoFiles)
        videoBase = strip_extension(videoFiles(j).name);
        videoNorm = normalize_name(videoBase);
        if strcmp(csvNorm, videoNorm)
            matchIdx = j;
            matchMethod = 'exact basename';
            break;
        end
    end

    if isempty(matchIdx)
        candidates = [];
        for j = 1:numel(videoFiles)
            videoBase = strip_extension(videoFiles(j).name);
            videoNorm = normalize_name(videoBase);
            if contains(csvNorm, videoNorm) || contains(videoNorm, csvNorm)
                candidates(end+1) = j; %#ok<AGROW>
            end
        end
        if isscalar(candidates)
            matchIdx = candidates;
            matchMethod = 'shared basename';
        end
    end

    if isempty(matchIdx) && numel(csvFiles) == numel(videoFiles)
        matchIdx = i;
        matchMethod = 'index fallback';
        warning('Using index-based match for %s -> %s.', csvFiles(i).name, videoFiles(i).name);
    end

    if isempty(matchIdx)
        warning('No MP4 match found for CSV file: %s', csvFiles(i).name);
        continue;
    end

    if usedVideos(matchIdx)
        warning('Video file was matched more than once: %s', videoFiles(matchIdx).name);
    end
    usedVideos(matchIdx) = true;

    pairs(end+1).csvFile = csvFiles(i).name; %#ok<AGROW>
    pairs(end).videoFile = videoFiles(matchIdx).name;
    pairs(end).csvPath = fullfile(csvFiles(i).folder, csvFiles(i).name);
    pairs(end).videoPath = fullfile(videoFiles(matchIdx).folder, videoFiles(matchIdx).name);
    pairs(end).baseName = strip_extension(videoFiles(matchIdx).name);
    pairs(end).matchMethod = matchMethod;
end
end

function out = strip_extension(fileName)
[~, out, ~] = fileparts(fileName);
end

function out = normalize_name(name)
out = lower(char(name));
out = regexprep(out, 'dlc.*$', '');
out = regexprep(out, '[^a-z0-9]', '');
end
