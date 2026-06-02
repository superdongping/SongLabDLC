clear all; 
close all;
clc;


% Search for XLSX files in the current folder
xlsx_files = dir('*.xlsx');
xlsx_files_name = xlsx_files(1).name;
fprintf('Reading xlsx file: %s\n', xlsx_files_name);

% Define the sheet number
sheet = 1;

% Read the data from the specified sheet
Raw_Data = readtable(xlsx_files_name, 'Sheet', sheet);

% Extract DLC data and manual data (assuming these are in columns 12 and 13)
dlc_data = Raw_Data{:, 17};  % Column 12: DLC data
manual_data = Raw_Data{:, 20};  % Column 13: Manual data

% Remove NaN values if necessary (optional, depending on your data)
valid_idx = ~isnan(dlc_data) & ~isnan(manual_data);
dlc_data = dlc_data(valid_idx);
manual_data = manual_data(valid_idx);

% Calculate Pearson correlation
correlation = corr(dlc_data, manual_data);
disp(['Pearson correlation: ', num2str(correlation)]);

% Fisher z-transformation for confidence interval calculation
n = length(dlc_data);  % Sample size
z = 0.5 * log((1 + correlation) / (1 - correlation));  % Fisher z-transformation
se = 1 / sqrt(n - 3);  % Standard error
z_conf_interval = [z - 1.96 * se, z + 1.96 * se];  % 95% confidence interval

% Convert z-values back to correlation scale
conf_interval = (exp(2 * z_conf_interval) - 1) ./ (exp(2 * z_conf_interval) + 1);
disp(['95% confidence interval for Pearson correlation: ', num2str(conf_interval(1)), ' to ', num2str(conf_interval(2))]);

% Plot Pearson correlation (scatter plot)
figure;
scatter(dlc_data, manual_data, 'filled');
xlabel('DLC Data');
ylabel('Manual Data');
title('Scatter Plot of DLC Data vs Manual Data with 95% CI');
hold on;

% Add the line of best fit
p = polyfit(dlc_data, manual_data, 1);
yfit = polyval(p, dlc_data);
plot(dlc_data, yfit, '-r', 'LineWidth', 1);

% Add 95% confidence interval band (optional, approximate visualization)
% Simulating the confidence intervals as a band around the fitted line
x_fit = sort(dlc_data);
y_fit_upper = polyval(p, x_fit) + (conf_interval(2) - correlation) * std(manual_data);  % Upper bound
y_fit_lower = polyval(p, x_fit) + (conf_interval(1) - correlation) * std(manual_data);  % Lower bound

% Plot the confidence interval band
fill([x_fit; flipud(x_fit)], [y_fit_lower; flipud(y_fit_upper)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Display the Pearson correlation and confidence interval on the plot
text(min(dlc_data), max(manual_data), ['r = ', num2str(correlation)], 'FontSize', 12, 'Color', 'blue');
text(min(dlc_data), max(manual_data) - 0.1 * (max(manual_data) - min(manual_data)), ...
    ['95% CI: [', num2str(conf_interval(1)), ', ', num2str(conf_interval(2)), ']'], 'FontSize', 12, 'Color', 'blue');

hold off;

