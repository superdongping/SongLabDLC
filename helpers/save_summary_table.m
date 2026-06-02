function outputPath = save_summary_table(T, outputDir, fileName, varargin)
%SAVE_SUMMARY_TABLE Save a summary table in the run output folder.

outputPath = fullfile(outputDir, fileName);
writetable(T, outputPath, varargin{:});
end
