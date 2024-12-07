% Define the paths
input_file_path_template = '/Users/lichaohui/Desktop/calculation/construction/r_forecast/regions/%s/regression_%s.xlsx';
output_figures_path_template = '/Users/lichaohui/Desktop/calculation/figures/%s';
output_combined_figure_path = '/Users/lichaohui/Desktop/calculation/construction/figures/combined_carbon_footprint.svg';

% Define the variables to analyze
variables = {'carbon'};
variable_titles = {'Carbon Footprint'};

% Define the countries to analyze
countries = {'China', 'France', 'Germany', 'India', 'Japan', 'US'};

% Define custom color schemes for each scenario
base_colors = {
    [0.8500, 0.3250, 0.0980]; % Base color for first scenario
    [0.3010, 0.7450, 0.9330]; % Base color for second scenario
    [0.4660, 0.6740, 0.1880]; % Base color for third scenario
    [0.9290, 0.6940, 0.1250]; % Base color for fourth scenario
    [0.4940, 0.1840, 0.5560]; % Base color for fifth scenario
};
%    [0.4940, 0.1840, 0.5560]; % Base color for fifth scenario

% Create a new figure for the combined plot with a larger size
figure('Position', [100, 100, 1400, 900]);

% Plot each country's carbon footprint in its designated subplot
for i = 1:length(countries)
    subplot(3, 2, i);
    input_file_path = sprintf(input_file_path_template, countries{i}, countries{i});
    plotCarbonFootprint(input_file_path, 'carbon', 'Carbon Footprint', base_colors, countries{i});
    title(countries{i});
end

% Save the combined figure
saveas(gcf, output_combined_figure_path);

% Function to plot carbon footprint for a given data path
function plotCarbonFootprint(data_path, var_name, var_title, base_colors, country)
    % Read the data from the Excel file
    opts = detectImportOptions(data_path, 'Sheet', 'merge');
    regression_data = readtable(data_path, opts);

    % Fill missing values with linear interpolation
    pop_columns = {'pop', 'pop_ssp1', 'pop_ssp2', 'pop_ssp3', 'pop_ssp4', 'pop_ssp5'};
    for col = 1:length(pop_columns)
        regression_data.(pop_columns{col}) = fillmissing(regression_data.(pop_columns{col}), 'linear');
    end

    % Extract variables for linear regression
    year = regression_data.year;
    target_var = regression_data.(var_name);
    pop = regression_data.pop;
    pop_ssp1 = regression_data.pop_ssp1;
    pop_ssp2 = regression_data.pop_ssp2;
    pop_ssp3 = regression_data.pop_ssp3;
    pop_ssp4 = regression_data.pop_ssp4;
    pop_ssp5 = regression_data.pop_ssp5;

    % Perform linear regression
    model1 = fitlm(pop, target_var);
    model_ssp1 = fitlm(pop_ssp1, target_var);
    model_ssp2 = fitlm(pop_ssp2, target_var);
    model_ssp3 = fitlm(pop_ssp3, target_var);
    model_ssp4 = fitlm(pop_ssp4, target_var);
    model_ssp5 = fitlm(pop_ssp5, target_var);

    % Predict future values
    predicted_var = predict(model1, pop);
    predicted_var_ssp1 = predict(model_ssp1, pop_ssp1);
    predicted_var_ssp2 = predict(model_ssp2, pop_ssp2);
    predicted_var_ssp3 = predict(model_ssp3, pop_ssp3);
    predicted_var_ssp4 = predict(model_ssp4, pop_ssp4);
    predicted_var_ssp5 = predict(model_ssp5, pop_ssp5);

    % Calculate 95% prediction intervals
    [pred_var, CI_var] = predict(model1, pop, 'Prediction', 'observation');
    [pred_var_ssp1, CI_var_ssp1] = predict(model_ssp1, pop_ssp1, 'Prediction', 'observation');
    [pred_var_ssp2, CI_var_ssp2] = predict(model_ssp2, pop_ssp2, 'Prediction', 'observation');
    [pred_var_ssp3, CI_var_ssp3] = predict(model_ssp3, pop_ssp3, 'Prediction', 'observation');
    [pred_var_ssp4, CI_var_ssp4] = predict(model_ssp4, pop_ssp4, 'Prediction', 'observation');
    [pred_var_ssp5, CI_var_ssp5] = predict(model_ssp5, pop_ssp5, 'Prediction', 'observation');

    % Define shades for each color (darker for central, lighter for bounds)
    shade_factor = 0.7;
    lighter_shade = @(color) color + (1 - color) * shade_factor;

    % Plotting the results
    hold on;

    % Plot prediction intervals with shades
    fill([year; flipud(year)], [CI_var(:,1); flipud(CI_var(:,2))], lighter_shade(base_colors{1}), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([year; flipud(year)], [CI_var_ssp1(:,1); flipud(CI_var_ssp1(:,2))], lighter_shade(base_colors{1}), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([year; flipud(year)], [CI_var_ssp2(:,1); flipud(CI_var_ssp2(:,2))], lighter_shade(base_colors{2}), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([year; flipud(year)], [CI_var_ssp3(:,1); flipud(CI_var_ssp3(:,2))], lighter_shade(base_colors{3}), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([year; flipud(year)], [CI_var_ssp4(:,1); flipud(CI_var_ssp4(:,2))], lighter_shade(base_colors{4}), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([year; flipud(year)], [CI_var_ssp5(:,1); flipud(CI_var_ssp5(:,2))], lighter_shade(base_colors{5}), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Plot fitted historical values (1995-2022) as pink lines
    h1 = plot(year(1:find(year == 2022)), predicted_var(1:find(year == 2022)), 'Color', [1, 0.6, 0.8], 'LineWidth', 2);

    % Plot projected values (2022 onwards) as median line
    h2 = plot(year(find(year == 2022):end), predicted_var(find(year == 2022):end), 'Color', base_colors{1}, 'LineWidth', 2);
    h3 = plot(year(find(year == 2022):end), predicted_var_ssp1(find(year == 2022):end), 'Color', base_colors{1}, 'LineWidth', 2);
    h4 = plot(year(find(year == 2022):end), predicted_var_ssp2(find(year == 2022):end), 'Color', base_colors{2}, 'LineWidth', 2);
    h5 = plot(year(find(year == 2022):end), predicted_var_ssp3(find(year == 2022):end), 'Color', base_colors{3}, 'LineWidth', 2);
    h6 = plot(year(find(year == 2022):end), predicted_var_ssp4(find(year == 2022):end), 'Color', base_colors{4}, 'LineWidth', 2);
    h7 = plot(year(find(year == 2022):end), predicted_var_ssp5(find(year == 2022):end), 'Color', base_colors{5}, 'LineWidth', 2);

    % Plot actual values as black lines for historical data
    h8 = plot(year(1:find(year == 2022)), target_var(1:find(year == 2022)), 'k-', 'LineWidth', 1.5);

    % Add labels and title with larger fonts
    title(var_title, 'FontSize', 14);

    % Customize plot
    set(gcf, 'Color', 'w');
    xlim([1995, 2050]);

    % Remove the grid and add a box around the plot
    box on;
    grid off;

    % Set custom ylim based on the country
    if strcmp(country, 'US')
        ylim([0, 14e11]);
    elseif strcmp(country, 'Germany')
        ylim([0, 4e11]);
    elseif strcmp(country, 'France')
        ylim([0, 3e11]);
    end
    
    % Adjust legend to reflect correct line colors
    legend([h8, h1, h2, h3, h4, h5, h6, h7], {'Historical', 'Fitted Historical', 'Predicted Median', 'Predicted SSP1', ...
        'Predicted SSP2', 'Predicted SSP3', 'Predicted SSP4', 'Predicted SSP5'}, ...
        'Location', 'eastoutside', 'FontSize', 10);
    
    hold off;
end
