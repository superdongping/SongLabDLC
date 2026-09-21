function video = open_behavior_video(path)
%OPEN_BEHAVIOR_VIDEO Native VideoReader or verified no-reencoding fallback.
readablePath = prepare_behavior_video(path);
video = VideoReader(readablePath);
end
