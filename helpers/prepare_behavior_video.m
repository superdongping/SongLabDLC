function [pathOut, method] = prepare_behavior_video(videoPath, forceRemux)
%PREPARE_BEHAVIOR_VIDEO Preflight native reader, then verified stream-copy MP4.
% Set SONGLABDLC_FFMPEG to ffmpeg.exe if it is not on PATH.
% No re-encoding, frame-rate conversion, resizing or source replacement.
if nargin < 2, forceRemux=false; end
persistent keys paths
if isempty(keys), keys={}; paths={}; end
videoPath=char(videoPath); info=dir(videoPath);
if isempty(info), error('SongLabDLC:MissingVideo','Video not found: %s',videoPath); end
key=sprintf('%s|%d|%.12f',videoPath,info.bytes,info.datenum);
idx=find(strcmp(keys,key),1);
if ~isempty(idx) && isfile(paths{idx})
    pathOut=paths{idx}; method='verified cached remux'; return;
end
try
    if forceRemux, error('Forced remux verification test.'); end
    v=VideoReader(videoPath); %#ok<TNMLP>
    if ~hasFrame(v), error('No video frames.'); end
    readFrame(v); clear v
    pathOut=videoPath; method='native'; return;
catch nativeError
    clear v
end
ffmpeg=getenv('SONGLABDLC_FFMPEG'); if isempty(ffmpeg), ffmpeg='ffmpeg'; end
folder=tempname; mkdir(folder);
pathOut=fullfile(folder,'compatible.mp4');
cleanup=onCleanup(@() cleanup_failed(folder));
cmd=sprintf('%s -hide_banner -v error -i %s -map 0:v:0 -an -c:v copy -movflags +faststart -n %s', ...
    quoted(ffmpeg),quoted(videoPath),quoted(pathOut));
[code,msg]=system(cmd);
if code~=0
    error('SongLabDLC:VideoDecode','Native read failed: %s\nStream-copy fallback failed: %s\nInstall FFmpeg or set SONGLABDLC_FFMPEG. Source retained.',nativeError.message,msg);
end
sourceHash=fullfile(folder,'source.txt'); destHash=fullfile(folder,'copy.txt');
framehash(ffmpeg,videoPath,sourceHash); framehash(ffmpeg,pathOut,destHash);
a=readhash(sourceHash); b=readhash(destHash);
if height(a)~=height(b) || isempty(a) || ~isequal(a.hash,b.hash) ...
        || any(abs((a.pts-a.pts(1))-(b.pts-b.pts(1)))>1100)
    error('SongLabDLC:RemuxMismatch','Frame content/count or relative timestamps changed. Original retained.');
end
v=VideoReader(pathOut); readFrame(v); clear v
% Mark successful cache before cleanup runs; source is never touched.
fid=fopen(fullfile(folder,'verified'),'w'); fclose(fid);
keys{end+1}=key; paths{end+1}=pathOut;
method='verified stream-copy MP4';
warning('SongLabDLC:Remux','Using verified temporary MP4 for %s. Original retained.',videoPath);
end
function q=quoted(s)
s=char(s);
if any(ismember(s,[char(34) char(10) char(13) '%!&|<>^']))
    error('SongLabDLC:UnsafePath','Path contains unsupported command-shell characters.');
end
q=[char(34) s char(34)];
end
function framehash(ffmpeg,source,destination)
cmd=sprintf('%s -hide_banner -v error -xerror -i %s -map 0:v:0 -an -fps_mode passthrough -enc_time_base 1:1000000 -f framemd5 -n %s',quoted(ffmpeg),quoted(source),quoted(destination));
[code,msg]=system(cmd);
if code~=0, error('SongLabDLC:VerificationFailed','Frame verification failed: %s',msg); end
end
function t=readhash(path)
lines=splitlines(string(fileread(path))); lines=lines(strlength(lines)>0 & ~startsWith(lines,'#'));
parts=split(lines,',');
t=table(str2double(strtrim(parts(:,3))),strtrim(parts(:,6)),'VariableNames',{'pts','hash'});
end
function cleanup_failed(folder)
if ~isfile(fullfile(folder,'verified')) && isfolder(folder)
    rmdir(folder,'s'); % Only this function's fresh tempname directory.
end
end
