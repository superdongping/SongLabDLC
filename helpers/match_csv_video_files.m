function pairs = match_csv_video_files(dataFolder, preferredExtension)
%MATCH_CSV_VIDEO_FILES Match DLC files by basename; never by directory order.
% preferredExtension (optional): e.g. '.mp4' resolves same-basename copies.
if nargin < 2, preferredExtension = ''; end
supported = {'.mp4','.avi','.mkv','.mov','.m4v','.mpg','.mpeg','.wmv','.webm','.mj2'};
files = dir(dataFolder); files = files(~[files.isdir]);
isVideo = false(size(files)); isCSV = false(size(files));
for i=1:numel(files)
    [~,base,ext] = fileparts(files(i).name);
    isVideo(i) = any(strcmpi(ext,supported)) && ~endsWith(base,'.pending','IgnoreCase',true);
    isCSV(i) = strcmpi(ext,'.csv') && ~startsWith(base,'summary_','IgnoreCase',true) ...
        && ~endsWith(base,'_frame_timing','IgnoreCase',true) && ~strcmpi(base,'SongScope_record_log');
end
videos = files(isVideo); csvs = files(isCSV);
if isempty(csvs), error('SongLabDLC:NoCSV','No candidate DLC CSV files in %s.',dataFolder); end
if isempty(videos), error('SongLabDLC:NoVideo','No supported video files in %s.',dataFolder); end
pairs = struct('csvFile',{},'videoFile',{},'csvPath',{},'videoPath',{},'baseName',{},'matchMethod',{});
used = false(size(videos));
for i=1:numel(csvs)
    [~,base] = fileparts(csvs(i).name);
    % DLC appends its model descriptor to the exact source-video basename.
    base = regexprep(base,'DLC.*$','','ignorecase');
    candidates = [];
    for j=1:numel(videos)
        [~,vbase] = fileparts(videos(j).name);
        if strcmpi(base,vbase), candidates(end+1)=j; end %#ok<AGROW>
    end
    if numel(candidates)>1 && ~isempty(preferredExtension)
        keep=false(size(candidates));
        for k=1:numel(candidates)
            [~,~,ext]=fileparts(videos(candidates(k)).name);
            keep(k)=strcmpi(ext,preferredExtension);
        end
        candidates=candidates(keep);
    end
    if numel(candidates)>1
        error('SongLabDLC:AmbiguousVideo', ...
            'Multiple videos match %s. Specify preferredExtension (e.g. ''.mp4'') or select a folder with one format per basename.',csvs(i).name);
    elseif isempty(candidates)
        warning('SongLabDLC:UnmatchedCSV','No exact video basename match for %s; skipped.',csvs(i).name);
        continue;
    end
    j=candidates(1);
    if used(j)
        error('SongLabDLC:DuplicateCSV','Multiple DLC CSV files match %s. Keep one chosen tracking result per video.',videos(j).name);
    end
    used(j)=true;
    [~,vbase]=fileparts(videos(j).name);
    pairs(end+1)=struct('csvFile',csvs(i).name,'videoFile',videos(j).name, ...
        'csvPath',fullfile(dataFolder,csvs(i).name),'videoPath',fullfile(dataFolder,videos(j).name), ...
        'baseName',vbase,'matchMethod','exact basename (DLC suffix removed)'); %#ok<AGROW>
end
end
