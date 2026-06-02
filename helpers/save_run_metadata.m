function metadataPath = save_run_metadata(runInfo, outputDir)
%SAVE_RUN_METADATA Save run metadata as JSON.

metadataPath = fullfile(outputDir, 'run_metadata.json');
txt = jsonencode(runInfo, 'PrettyPrint', true);

fid = fopen(metadataPath, 'w');
if fid == -1
    error('Could not create metadata file: %s', metadataPath);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
clear cleanupObj;
end
