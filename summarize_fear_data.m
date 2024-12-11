% summarize_fear_data.m
% This script parses multiple .txt files containing fear condition test results
% and compiles selected data into a summary Excel file.

% Clear workspace and command window for cleanliness
clear;
clc;

% Define the name of the summary Excel file
summaryFileName = 'summary_fear.xlsx';

% Get a list of all .txt files in the current directory
txtFiles = dir('*.txt');

% Initialize arrays to store extracted data
experimentIDs      = {}; % Changed to cell array to handle strings
totalTimes         = [];
avgMotions         = [];
freezeEpisodes     = [];
timeFreezes        = [];
pctFreezings       = [];

% Loop through each .txt file
for i = 1:length(txtFiles)
    % Get the current file name
    currentFile = txtFiles(i).name;
    
    % Open the file for reading
    fid = fopen(currentFile, 'r');
    if fid == -1
        warning('Could not open file: %s. Skipping...', currentFile);
        continue; % Skip to the next file if unable to open
    end
    
    % Read all lines from the file
    fileLines = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid); % Close the file after reading
    fileLines = fileLines{1}; % Extract the cell array of lines
    
    % Initialize variables for the current file
    experimentID      = ''; % Changed to string
    totalTime         = NaN;
    avgMotion         = NaN;
    freezeEpisode     = NaN;
    timeFreeze        = NaN;
    pctFreezing       = NaN;
    
    % === Extract Experiment ID ===
    % Look for the line that contains 'Experiment'
    expLineIdx = find(contains(fileLines, 'Experiment'), 1);
    if ~isempty(expLineIdx)
        % Use regular expression to extract the text after ':'
        % This handles both numerical and alphanumeric IDs
        tokens = regexp(fileLines{expLineIdx}, 'Experiment\s*:\s*(\S+)', 'tokens');
        if ~isempty(tokens)
            experimentID = tokens{1}{1}; % Capture as string
        else
            warning('Experiment ID not found or invalid in file: %s', currentFile);
            continue; % Skip if Experiment ID is invalid
        end
    else
        warning('Experiment line not found in file: %s', currentFile);
        continue; % Skip if Experiment line is missing
    end
    
    % === Extract Metrics ===
    % Look for the line that contains 'Total Time (s)'
    metricsHeaderIdx = find(contains(fileLines, 'Total Time (s)'), 1);
    if ~isempty(metricsHeaderIdx)
        % The next line after the header contains the data
        dataLineIdx = metricsHeaderIdx + 2; % Assuming there's a divider line
        if dataLineIdx <= length(fileLines)
            dataLine = strtrim(fileLines{dataLineIdx});
            % Split the data line by whitespace or tabs
            dataParts = regexp(dataLine, '\s+', 'split');
            dataParts = dataParts(~cellfun('isempty', dataParts)); % Remove empty cells
            
            if length(dataParts) >= 5
                % Convert each part to a number
                totalTime       = str2double(dataParts{1});
                avgMotion       = str2double(dataParts{2});
                freezeEpisode   = str2double(dataParts{3});
                timeFreeze      = str2double(dataParts{4});
                pctFreezing     = str2double(dataParts{5});
            else
                warning('Insufficient data fields in file: %s', currentFile);
                continue; % Skip if data fields are insufficient
            end
        else
            warning('Data line missing after metrics header in file: %s', currentFile);
            continue; % Skip if data line is missing
        end
    else
        warning('Metrics header not found in file: %s', currentFile);
        continue; % Skip if metrics header is missing
    end
    
    % === Append Extracted Data ===
    experimentIDs{end+1,1}    = experimentID;      %#ok<*AGROW>
    totalTimes(end+1,1)       = totalTime;
    avgMotions(end+1,1)       = avgMotion;
    freezeEpisodes(end+1,1)   = freezeEpisode;
    timeFreezes(end+1,1)      = timeFreeze;
    pctFreezings(end+1,1)     = pctFreezing;
end

% === Create a Table from the Extracted Data ===
summaryTable = table(...
    experimentIDs, ...
    totalTimes, ...
    avgMotions, ...
    freezeEpisodes, ...
    timeFreezes, ...
    pctFreezings, ...
    'VariableNames', {
        'Experiment_ID', ...
        'Total_Time_s', ...
        'Avg_Motion', ...
        'Freeze_Episodes', ...
        'Time_Freeze_s', ...
        'Pct_Freezing' ...
    });

% === Write the Table to an Excel File ===
writetable(summaryTable, summaryFileName);

% Display a message upon successful completion
fprintf('Summary data has been successfully written to %s\n', summaryFileName);
