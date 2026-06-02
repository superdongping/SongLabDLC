function idx = limit_tracking_duration(frames, fps, durationSec)
%LIMIT_TRACKING_DURATION Logical index for a maximum analysis duration.

if isnan(durationSec) || isempty(durationSec)
    idx = true(size(frames));
    return;
end

maxFrame = round(durationSec * fps);
idx = frames <= maxFrame;
end
