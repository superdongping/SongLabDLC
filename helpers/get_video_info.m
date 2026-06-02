function info = get_video_info(videoPath)
%GET_VIDEO_INFO Read basic metadata from an MP4 video.

videoObj = VideoReader(videoPath);
info = struct();
info.path = videoPath;
info.file = string(get_file_name(videoPath));
info.frame_rate = videoObj.FrameRate;
info.duration_sec = videoObj.Duration;
info.width = videoObj.Width;
info.height = videoObj.Height;
info.estimated_num_frames = floor(videoObj.Duration * videoObj.FrameRate);
end

function name = get_file_name(path)
[~, base, ext] = fileparts(path);
name = [base ext];
end
