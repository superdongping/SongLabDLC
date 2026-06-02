function options = get_default_behavior_options(assayKey)
%GET_DEFAULT_BEHAVIOR_OPTIONS Centralized assay parameters.

options = struct();
options.assay = assayKey;
options.analysis_duration_sec = NaN;
options.likelihood_threshold = 0.05;
options.body_part_main = {'MainBody', 'Mainbody', 'Body', 'Center'};
options.body_part_nose = {'Nose'};
options.body_part_head = {'Head'};

switch assayKey
    case 'OFT'
        options.analysis_duration_sec = 6 * 60;
        options.likelihood_threshold = 0.05;
        options.open_field_scale_m = 0.45;
        options.center_square_m = 0.25;
        options.heatmap_grid_size = [100, 100];
        options.heatmap_sigma = 2;
        options.heatmap_display_max_fraction = 0.20;
    case 'NPR'
        options.analysis_duration_sec = 6 * 60;
        options.likelihood_threshold = 0.03;
        options.open_field_scale_m = 0.45;
        options.object_roi_radius_m = 0.04;
        options.window_duration_sec = 0.2;
        options.min_overlap = 0.8;
        options.movement_threshold_px = 1.5;
        options.sitting_stability_sec = 2;
        options.merge_gap_sec = 0.40;
        options.min_bout_sec = 0.15;
    case 'ZERO_MAZE'
        options.analysis_duration_sec = 6 * 60;
        options.likelihood_threshold = 0.01;
        options.heatmap_grid_bins = 100;
        options.heatmap_sigma = 2;
    case 'Y_MAZE'
        options.analysis_duration_sec = 8 * 60;
        options.likelihood_threshold = 0.01;
    case 'FST'
        options.analysis_duration_sec = 300;
        options.scale_height_cm = 24.25;
        options.speed_threshold_cm_s = 20;
        options.continuous_duration_sec = 1;
        options.immobility_threshold_percent = 70;
        options.body_parts = {'Head', 'MainBody', 'Butt', 'TailTip', ...
            'LeftForehand', 'RightForehand', 'LeftHindpaw', 'RightHindpaw'};
        options.weights = struct('Head', 5, 'MainBody', 5, 'Butt', 5, 'TailTip', 5, ...
            'LeftForehand', 10, 'RightForehand', 10, 'LeftHindpaw', 30, 'RightHindpaw', 30);
    case 'TST'
        options.analysis_duration_sec = 360;
        options.scale_height_cm = 24.25;
        options.speed_threshold_cm_s = 10;
        options.continuous_duration_sec = 1;
        options.immobility_threshold_percent = 70;
        options.body_parts = {'Head', 'MainBody', 'Butt', 'TailTip', ...
            'LeftForehand', 'RightForehand', 'LeftHindpaw', 'RightHindpaw'};
        options.weights = struct('Head', 30, 'MainBody', 20, 'Butt', 10, 'TailTip', 0, ...
            'LeftForehand', 10, 'RightForehand', 10, 'LeftHindpaw', 10, 'RightHindpaw', 10);
    otherwise
        error('No default options defined for assay: %s', assayKey);
end
end
