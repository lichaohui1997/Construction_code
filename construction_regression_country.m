% Define the paths
input_file_path_template = '/Users/lichaohui/Desktop/calculation/construction/r_forecast/regions/%s/regression_%s.xlsx';
output_file_path = '/Users/lichaohui/Desktop/calculation/construction/SIdata3.xlsx';
output_figures_path_template = '/Users/lichaohui/Desktop/calculation/construction/figures/%s';

% Define the variables to analyze
variables = {'carbon', 'energy', 'water'};
variable_titles = {'Carbon Footprint', 'Energy Footprint', 'Water Footprint'};

% Define the countries to analyze
countries = {'China', 'France', 'Germany', 'India', 'Japan', 'US'};

% Loop through each country
for country_idx = 1:length(countries)
    country = countries{country_idx};
    
    % Loop through each variable
    for var_idx = 1:length(variables)
        var_name = variables{var_idx};
        var_title = variable_titles{var_idx};
        
        % Generate input file path
        input_file_path = sprintf(input_file_path_template, country, country);
        
        % Read the data from the Excel file
        opts = detectImportOptions(input_file_path, 'Sheet', 'merge');
        regression_data = readtable(input_file_path, opts);
        
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
        
        % Define custom colors
        color1 = [0, 0.4470, 0.7410]; % Blue
        color2 = [0.8500, 0.3250, 0.0980]; % Red
        color3 = [0.9290, 0.6940, 0.1250]; % Yellow
        color4 = [0.4940, 0.1840, 0.5560]; % Purple
        color5 = [0.4660, 0.6740, 0.1880]; % Green
        color6 = [0.3010, 0.7450, 0.9330]; % Cyan
        
        % Plotting the results
        figure;
        hold on;
        
        % Plot prediction intervals
        fill([year; flipud(year)], [CI_var(:,1); flipud(CI_var(:,2))], color1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill([year; flipud(year)], [CI_var_ssp1(:,1); flipud(CI_var_ssp1(:,2))], color2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill([year; flipud(year)], [CI_var_ssp2(:,1); flipud(CI_var_ssp2(:,2))], color3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill([year; flipud(year)], [CI_var_ssp3(:,1); flipud(CI_var_ssp3(:,2))], color4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill([year; flipud(year)], [CI_var_ssp4(:,1); flipud(CI_var_ssp4(:,2))], color5, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        fill([year; flipud(year)], [CI_var_ssp5(:,1); flipud(CI_var_ssp5(:,2))], color6, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Plot predicted values
        plot(year, predicted_var, 'Color', color1, 'LineWidth', 1.5);
        plot(year, predicted_var_ssp1, 'Color', color2, 'LineWidth', 1.5);
        plot(year, predicted_var_ssp2, 'Color', color3, 'LineWidth', 1.5);
        plot(year, predicted_var_ssp3, 'Color', color4, 'LineWidth', 1.5);
        plot(year, predicted_var_ssp4, 'Color', color5, 'LineWidth', 1.5);
        plot(year, predicted_var_ssp5, 'Color', color6, 'LineWidth', 1.5);
        
        % Plot actual values
        plot(year, target_var, 'k.', 'MarkerSize', 10);
        
        % Add labels and title
        title(['Predicted ', var_title, ' of Construction Sector in ', country]);
        xlabel('Year');
        ylabel([var_title]);
        
        % Customize plot
        set(gcf, 'Color', 'w');
        grid on;
        legend({'Prediction Interval 1', 'Prediction Interval SSP1', 'Prediction Interval SSP2', ...
            'Prediction Interval SSP3', 'Prediction Interval SSP4', 'Prediction Interval SSP5', ...
            'Predicted', 'Predicted SSP1', 'Predicted SSP2', ...
            'Predicted SSP3', 'Predicted SSP4', 'Predicted SSP5', 'Actual'}, ...
            'Location', 'Best');
        
        hold off;
        
        % Generate output figure path
        output_figures_path = sprintf(output_figures_path_template, country);
        if ~exist(output_figures_path, 'dir')
            mkdir(output_figures_path);
        end
        
        % Save the plot
        saveas(gcf, fullfile(output_figures_path, ['predicted_', var_name, '_footprint_ssps.svg']));
        
        % Export the variable "regression_data" to the output Excel file
        writetable(regression_data, output_file_path, 'Sheet', [country, '_', var_name]);
    end
end
