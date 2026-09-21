function SongLabDLC_behavior_analysis(assayName, dataFolder, preferredExtension)
%SONGLABDLC_BEHAVIOR_ANALYSIS Unified entry point for post-DLC behavior analysis.
%
% Run with no inputs for an interactive workflow:
%   SongLabDLC_behavior_analysis
%
% Or provide inputs programmatically:
%   SongLabDLC_behavior_analysis("OFT", "C:\path\to\data", ".mp4")

interactive = nargin < 2;
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
addpath(fullfile(rootDir, 'helpers'));
addpath(fullfile(rootDir, 'assays'));

assays = ["OFT", "NPR", "Zero Maze", "Y-maze", "FST", "TST"];

if nargin < 1 || strlength(string(assayName)) == 0
    [idx, ok] = listdlg( ...
        'PromptString', 'Select behavior assay', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(assays), ...
        'Name', 'SongLabDLC Behavior Analysis');
    if ~ok
        disp('Analysis canceled.');
        return;
    end
    assayName = assays(idx);
else
    assayName = string(assayName);
end

assayKey = normalize_assay_name(assayName);

if nargin < 2 || strlength(string(dataFolder)) == 0
    dataFolder = uigetdir(pwd, 'Select folder containing DLC CSV and video files');
    if isequal(dataFolder, 0)
        disp('Analysis canceled.');
        return;
    end
end
dataFolder = char(dataFolder);

options = get_default_behavior_options(assayKey);
if nargin < 3, preferredExtension = ''; end
try
    pairs = match_csv_video_files(dataFolder, preferredExtension);
catch err
    if interactive && strcmp(err.identifier,'SongLabDLC:AmbiguousVideo')
        choices={'.mp4','.mkv','.avi','.mov','.m4v','.mpg','.mpeg','.wmv','.webm','.mj2'};
        [choice,ok]=listdlg('PromptString','Same-basename copies found. Choose the video format to analyze.', ...
            'SelectionMode','single','ListString',choices,'Name','Choose video format');
        if ~ok, disp('Analysis canceled.'); return; end
        preferredExtension=choices{choice};
        pairs = match_csv_video_files(dataFolder, preferredExtension);
    else
        rethrow(err);
    end
end
for i=1:numel(pairs)
    pairs(i).originalVideoPath = pairs(i).videoPath;
    [pairs(i).videoPath,pairs(i).readerMethod] = prepare_behavior_video(pairs(i).videoPath);
end

if isempty(pairs)
    error('No matched CSV/video pairs were found in: %s', dataFolder);
end

fprintf('\nMatched files:\n');
disp(struct2table(pairs));

outputRoot = fullfile(dataFolder, 'SongLabDLC_behavior_results');
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

runStamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
runOutputDir = fullfile(outputRoot, sprintf('%s_%s', assayKey, runStamp));
mkdir(runOutputDir);

runInfo = struct();
runInfo.pipeline = 'SongLabDLC behavior pipeline';
runInfo.pipeline_version = '0.2.0';
runInfo.preferred_video_extension = preferredExtension;
runInfo.timing_assumption = 'Existing assay metrics use nominal/average FPS, not per-frame timestamps. Review irregular-rate recordings before time-based interpretation.';
runInfo.assay = assayKey;
runInfo.data_folder = dataFolder;
runInfo.output_folder = runOutputDir;
runInfo.run_datetime = char(datetime('now'));
runInfo.options = options;
runInfo.file_pairs = pairs;

switch assayKey
    case 'OFT'
        result = analyze_oft(pairs, runOutputDir, options);
    case 'NPR'
        result = analyze_npr(pairs, runOutputDir, options);
    case 'ZERO_MAZE'
        result = analyze_zero_maze(pairs, runOutputDir, options);
    case 'Y_MAZE'
        result = analyze_y_maze(pairs, runOutputDir, options);
    case 'FST'
        result = analyze_fst(pairs, runOutputDir, options);
    case 'TST'
        result = analyze_tst(pairs, runOutputDir, options);
    otherwise
        error('Unsupported assay: %s', assayName);
end

runInfo.result = result;
metadataPath = save_run_metadata(runInfo, runOutputDir);

fprintf('\nAnalysis complete.\n');
fprintf('Output folder: %s\n', runOutputDir);
fprintf('Metadata: %s\n', metadataPath);
end

function key = normalize_assay_name(name)
name = upper(strtrim(string(name)));
name = replace(name, "-", "_");
name = replace(name, " ", "_");

switch name
    case {"OFT", "OPEN_FIELD", "OPEN_FIELD_TEST"}
        key = 'OFT';
    case {"NPR", "NOR", "NOVEL_OBJECT_RECOGNITION", "NOVELTY_PREFERENCE_RECOGNITION"}
        key = 'NPR';
    case {"ZERO_MAZE", "ZEROMAZE"}
        key = 'ZERO_MAZE';
    case {"Y_MAZE", "YMAZE"}
        key = 'Y_MAZE';
    case "FST"
        key = 'FST';
    case "TST"
        key = 'TST';
    otherwise
        error('Unknown assay name: %s', name);
end
end
