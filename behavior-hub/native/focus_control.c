/* Behavior Hub additions to FFmpeg 7.1.3 dshow.c, LGPL-2.1-or-later.
 * Control the SAME source filter that supplies captured packets. No second
 * camera instance, global driver assumption, or OpenCV side connection.
 */
static int bh_focus_setup(AVFormatContext *s)
{
    struct dshow_ctx *ctx = s->priv_data;
    IAMCameraControl *camera = NULL;
    long low, high, step, def, caps, value, flags, wanted;
    HRESULT hr;
    int ret = AVERROR(EIO);
    if (ctx->bh_focus_mode < 0)
        return 0;
    hr = IBaseFilter_QueryInterface(ctx->device_filter[VideoDevice],
                                    &IID_IAMCameraControl, (void **)&camera);
    if (FAILED(hr) || !camera) {
        av_log(s, AV_LOG_ERROR, "BH_FOCUS_ERROR Camera does not expose focus control.\n");
        return ret;
    }
    hr = IAMCameraControl_GetRange(camera, CameraControl_Focus,
                                  &low, &high, &step, &def, &caps);
    if (FAILED(hr) || step <= 0 || high < low) {
        av_log(s, AV_LOG_ERROR, "BH_FOCUS_ERROR Cannot read focus range.\n");
        goto end;
    }
    hr = IAMCameraControl_Get(camera, CameraControl_Focus, &value, &flags);
    if (FAILED(hr)) {
        av_log(s, AV_LOG_ERROR, "BH_FOCUS_ERROR Cannot read focus state.\n");
        goto end;
    }
    if (ctx->bh_focus_mode > 0) {
        wanted = ctx->bh_focus_mode == 1 ? CameraControl_Flags_Auto : CameraControl_Flags_Manual;
        if (!(caps & wanted)) {
            av_log(s, AV_LOG_ERROR, "BH_FOCUS_ERROR Requested focus mode is unsupported.\n");
            goto end;
        }
        if (ctx->bh_focus_mode == 2) {
            value = ctx->bh_focus_value;
            if (value < low || value > high || ((int64_t)value - low) % step) {
                av_log(s, AV_LOG_ERROR, "BH_FOCUS_ERROR Focus value is outside device range/step.\n");
                goto end;
            }
        }
        hr = IAMCameraControl_Set(camera, CameraControl_Focus, value, wanted);
        if (FAILED(hr)) {
            av_log(s, AV_LOG_ERROR, "BH_FOCUS_ERROR Device rejected focus setting.\n");
            goto end;
        }
        /* No frames pass the callback gate while the lens is settling. */
        Sleep(500);
        hr = IAMCameraControl_Get(camera, CameraControl_Focus, &value, &flags);
        if (FAILED(hr) || (flags & 3) != wanted ||
            (ctx->bh_focus_mode == 2 && value != ctx->bh_focus_value)) {
            av_log(s, AV_LOG_ERROR, "BH_FOCUS_ERROR Focus readback mismatch; capture blocked.\n");
            goto end;
        }
    }
    av_log(s, AV_LOG_INFO,
           "BH_FOCUS {\"minimum\":%ld,\"maximum\":%ld,\"step\":%ld,\"default\":%ld,"
           "\"capabilities\":%ld,\"value\":%ld,\"flags\":%ld,\"verified\":%s,"
           "\"source\":\"active_capture_filter\",\"protocol\":1}\n",
           low, high, step, def, caps, value, flags,
           ctx->bh_focus_mode > 0 ? "true" : "false");
    InterlockedExchange(&ctx->bh_focus_ready, 1);
    ret = 0;
end:
    IAMCameraControl_Release(camera);
    return ret;
}
