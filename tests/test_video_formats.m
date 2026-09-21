function test_video_formats(fixtureDir)
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'helpers'));
fprintf('MATLAB %s; %s\n',version,computer);
formats={'mp4','mkv','mov','avi'};
for i=1:numel(formats)
    file=fullfile(fixtureDir,['sample.' formats{i}]);
    v=open_behavior_video(file); n=0;
    while hasFrame(v), readFrame(v); n=n+1; end
    assert(n==25,'Unexpected decoded frame count for %s',formats{i});
    fprintf('READ_OK %s %d frames\n',formats{i},n); clear v
end
[path,method]=prepare_behavior_video(fullfile(fixtureDir,'sample.mkv'),true);
v=VideoReader(path); n=0; while hasFrame(v),readFrame(v);n=n+1;end
assert(n==25); fprintf('FALLBACK_OK %s\n',method);clear v
fid=fopen(fullfile(fixtureDir,'sampleDLC_model_filtered.csv'),'w');fprintf(fid,'fixture');fclose(fid);
try
    match_csv_video_files(fixtureDir); error('Expected conflict');
catch e
    assert(strcmp(e.identifier,'SongLabDLC:AmbiguousVideo'),e.message);
end
pairs=match_csv_video_files(fixtureDir,'.mp4');assert(numel(pairs)==1);assert(strcmp(pairs.videoFile,'sample.mp4'));
fid=fopen(fullfile(fixtureDir,'sample_frame_timing.csv'),'w');fprintf(fid,'not DLC');fclose(fid);
pairs=match_csv_video_files(fixtureDir,'.mkv');assert(numel(pairs)==1);
fid=fopen(fullfile(fixtureDir,'sampleDLC_other.csv'),'w');fprintf(fid,'fixture');fclose(fid);
try
    match_csv_video_files(fixtureDir,'.mp4');error('Expected duplicate');
catch e
    assert(strcmp(e.identifier,'SongLabDLC:DuplicateCSV'),e.message);
end
files=[dir(fullfile(root,'helpers','*.m'));dir(fullfile(root,'assays','*.m'))];
for i=1:numel(files)
    issues=checkcode(fullfile(files(i).folder,files(i).name),'-id');
    for j=1:numel(issues)
        assert(~strcmp(issues(j).id,'PARSE'),issues(j).message);
    end
end
fprintf('MATCH_AND_SYNTAX_OK\n');
end
