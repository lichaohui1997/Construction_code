% MATLAB Script to Create the Carbon Budget Pathways Plot with Historical Data

% Historical data from 1995 to 2022 (in GtCO2)
historical_years = 1995:2022;
historical_emissions = [29.2, 30.2, 31.9, 30.5, 30.8, 30.9, 30.7, 31.6, 33.5, ...
                        34.0, 34.5, 35.8, 36.1, 36.8, 36.7, 38.5, 39.7, 40.3, ...
                        40.1, 40.7, 41.1, 40.1, 40.6, 41.1, 41.6, 39.3, 41.1, 41.5];

% Projection settings
years_projection = 2022:1:2100;
y0_projection = 41.5;  % Start projection from 41.5 GtCO2 in 2023
RCBs = [2000, 1450, 1150, 950, 800];  % Remaining Carbon Budgets in GtCO2
probabilities = [17, 33, 50, 67, 83];  % Corresponding probabilities

% Initialize figure
figure;
hold on;
colors = {'r', 'g', 'b', 'm', 'c'};  % Color for each line

% Plot the historical data
plot(historical_years, historical_emissions, 'k', 'LineWidth', 2);
text(2010, 35, 'Historical', 'FontSize', 12, 'Color', 'k');

% Loop over each scenario and calculate the decay constant 'k'
for i = 1:length(RCBs)

    % Define the initial guess for k and the tolerance
    k_initial_guess = 0.1;
    tolerance = 1e-6;
    
    % Use fzero to find the correct k that results in the exact cumulative emissions (RCB)
    k = fzero(@(k) trapz(years_projection, y0_projection * exp(-k * (years_projection - 2022))) - RCBs(i), k_initial_guess);
    
    % Compute the path using the calculated k
    y_projection = y0_projection * exp(-k * (years_projection - 2022));

    % Plot the path
    plot(years_projection, y_projection, 'Color', colors{i}, 'LineWidth', 2);

    % Add labels to the lines
    text(2050, y_projection(find(years_projection==2050)), ...
        ['p=' num2str(probabilities(i)) '%, RCB=' num2str(RCBs(i)) 'GtCO2'], ...
        'FontSize', 12, 'Color', colors{i});
end

% Annotate the plot
xlabel('Year');
ylabel('GtCO_2');
title('Global CO_2 Pathways for 1.5°C using IPCC AR6 Remaining Carbon Budget');
legend('Historical Data', 'RCB=900GtCO2', 'RCB=650GtCO2', 'RCB=500GtCO2', ...
    'RCB=400GtCO2', 'RCB=300GtCO2', 'Location', 'northeast');
grid on;
xlim([1995 2050]);
ylim([0 45]);

hold off;

saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/carbon_budget.svg');

%% Combined Carbon Footprint and Carbon Budget Pathways Plot

% Load the Excel file
filename = '/Users/lichaohui/Desktop/calculation/construction/r_forecast/world/regression.xlsx';

% Read the data from the sheet named 'merge'
opts = detectImportOptions(filename, 'Sheet', 'merge');
regression_data = readtable(filename, opts);

% Extract variables for linear regression
year = regression_data.year;
carbon = regression_data.carbon / 1e12; % Convert to GtCO2
pop = regression_data.pop;
pop_ssp1 = regression_data.pop_ssp1;
pop_ssp2 = regression_data.pop_ssp2;
pop_ssp3 = regression_data.pop_ssp3;
pop_ssp4 = regression_data.pop_ssp4;
pop_ssp5 = regression_data.pop_ssp5;

% Perform linear regression
model1 = fitlm(pop, carbon);
model_ssp1 = fitlm(pop_ssp1, carbon);
model_ssp2 = fitlm(pop_ssp2, carbon);
model_ssp3 = fitlm(pop_ssp3, carbon);
model_ssp4 = fitlm(pop_ssp4, carbon);
model_ssp5 = fitlm(pop_ssp5, carbon);

% Predict future carbon emissions
predicted_carbon = predict(model1, pop);
predicted_carbon_ssp1 = predict(model_ssp1, pop_ssp1);
predicted_carbon_ssp2 = predict(model_ssp2, pop_ssp2);
predicted_carbon_ssp3 = predict(model_ssp3, pop_ssp3);
predicted_carbon_ssp4 = predict(model_ssp4, pop_ssp4);
predicted_carbon_ssp5 = predict(model_ssp5, pop_ssp5);

% Calculate 95% prediction intervals
[pred_carbon, CI_carbon] = predict(model1, pop, 'Prediction', 'observation');
[pred_carbon_ssp1, CI_carbon_ssp1] = predict(model_ssp1, pop_ssp1, 'Prediction', 'observation');
[pred_carbon_ssp2, CI_carbon_ssp2] = predict(model_ssp2, pop_ssp2, 'Prediction', 'observation');
[pred_carbon_ssp3, CI_carbon_ssp3] = predict(model_ssp3, pop_ssp3, 'Prediction', 'observation');
[pred_carbon_ssp4, CI_carbon_ssp4] = predict(model_ssp4, pop_ssp4, 'Prediction', 'observation');
[pred_carbon_ssp5, CI_carbon_ssp5] = predict(model_ssp5, pop_ssp5, 'Prediction', 'observation');

% Historical data from 1995 to 2022 (in GtCO2)
historical_years = 1995:2022;
historical_emissions = [29.2, 30.2, 31.9, 30.5, 30.8, 30.9, 30.7, 31.6, 33.5, ...
                        34.0, 34.5, 35.8, 36.1, 36.8, 36.7, 38.5, 39.7, 40.3, ...
                        40.1, 40.7, 41.1, 40.1, 40.6, 41.1, 41.6, 39.3, 41.1, 41.5];

% Projection settings for total carbon footprint
years_projection = 2022:1:2100;
y0_projection = 41.5;  % Start projection from 41.5 GtCO2 in 2022
RCBs = [2000, 1450, 1150, 950, 800];  % Remaining Carbon Budgets in GtCO2
probabilities = [17, 33, 50, 67, 83];  % Corresponding probabilities

% Initialize figure
figure;
hold on;

% Define custom colors
color1 = [0, 0.4470, 0.7410]; % Blue
color2 = [0.8500, 0.3250, 0.0980]; % Red
color3 = [0.9290, 0.6940, 0.1250]; % Yellow
color4 = [0.4940, 0.1840, 0.5560]; % Purple
color5 = [0.4660, 0.6740, 0.1880]; % Green
color6 = [0.3010, 0.7450, 0.9330]; % Cyan

% Plot prediction intervals for the construction carbon footprint
fill([year; flipud(year)], [CI_carbon(:,1); flipud(CI_carbon(:,2))], color1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([year; flipud(year)], [CI_carbon_ssp1(:,1); flipud(CI_carbon_ssp1(:,2))], color2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([year; flipud(year)], [CI_carbon_ssp2(:,1); flipud(CI_carbon_ssp2(:,2))], color3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([year; flipud(year)], [CI_carbon_ssp3(:,1); flipud(CI_carbon_ssp3(:,2))], color4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([year; flipud(year)], [CI_carbon_ssp4(:,1); flipud(CI_carbon_ssp4(:,2))], color5, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([year; flipud(year)], [CI_carbon_ssp5(:,1); flipud(CI_carbon_ssp5(:,2))], color6, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot predicted values for the construction carbon footprint
plot(year, predicted_carbon, 'Color', color1, 'LineWidth', 1.5);
plot(year, predicted_carbon_ssp1, 'Color', color2, 'LineWidth', 1.5);
plot(year, predicted_carbon_ssp2, 'Color', color3, 'LineWidth', 1.5);
plot(year, predicted_carbon_ssp3, 'Color', color4, 'LineWidth', 1.5);
plot(year, predicted_carbon_ssp4, 'Color', color5, 'LineWidth', 1.5);
plot(year, predicted_carbon_ssp5, 'Color', color6, 'LineWidth', 1.5);

% Plot actual carbon values for the construction carbon footprint
plot(year, carbon, 'k.', 'MarkerSize', 10);

% Plot the historical data for total carbon footprint
plot(historical_years, historical_emissions, 'k', 'LineWidth', 2);

% Loop over each scenario for the total carbon footprint
for i = 1:length(RCBs)

    % Define the initial guess for k and the tolerance
    k_initial_guess = 0.1;
    tolerance = 1e-6;
    
    % Use fzero to find the correct k that results in the exact cumulative emissions (RCB)
    k = fzero(@(k) trapz(years_projection, y0_projection * exp(-k * (years_projection - 2022))) - RCBs(i), k_initial_guess);
    
    % Compute the path using the calculated k
    y_projection = y0_projection * exp(-k * (years_projection - 2022));

    % Plot the path for total carbon footprint
    plot(years_projection, y_projection, 'Color', colors{i}, 'LineWidth', 2);

    % Add labels to the lines for total carbon footprint
    text(2050, y_projection(find(years_projection==2050)), ...
        ['p=' num2str(probabilities(i)) '%, RCB=' num2str(RCBs(i)) 'GtCO2'], ...
        'FontSize', 12, 'Color', colors{i});
end

% Add labels and title
title('Combined Predicted Carbon Footprint and Carbon Budget Pathways');
xlabel('Year');
ylabel('Carbon Footprint / Total CO_2 (GtCO_2)');

% Create a separate legend box
legend({'Construction Prediction Interval', 'Construction Prediction Interval SSP1', ...
    'Construction Prediction Interval SSP2', 'Construction Prediction Interval SSP3', ...
    'Construction Prediction Interval SSP4', 'Construction Prediction Interval SSP5', ...
    'Construction Predicted Carbon', 'Construction Predicted Carbon SSP1', ...
    'Construction Predicted Carbon SSP2', 'Construction Predicted Carbon SSP3', ...
    'Construction Predicted Carbon SSP4', 'Construction Predicted Carbon SSP5', ...
    'Actual Construction Carbon', 'Total Carbon Historical Data', 'Total Carbon RCB 500GtCO2', ...
    'Total Carbon RCB 300GtCO2', 'Total Carbon RCB 250GtCO2', ...
    'Total Carbon RCB 150GtCO2', 'Total Carbon RCB 100GtCO2'}, ...
    'Location', 'northeastoutside', 'Box', 'on');

% Customize plot
set(gcf, 'Color', 'w');
grid on;
xlim([1995 2050]);
ylim([0 45]);

hold off;

% Save the plot
saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/combined_carbon_footprint_and_budget.svg');
%% Final Simplified Combined Carbon Footprint and Carbon Budget Pathways Plot

% Load the Excel file
filename = '/Users/lichaohui/Desktop/calculation/construction/r_forecast/world/regression.xlsx';

% Read the data from the sheet named 'merge'
opts = detectImportOptions(filename, 'Sheet', 'merge');
regression_data = readtable(filename, opts);

% Extract variables for linear regression
year = regression_data.year;
carbon = regression_data.carbon / 1e12; % Convert to GtCO2
pop = regression_data.pop;
pop_ssp2 = regression_data.pop_ssp2;

% Perform linear regression
model_ssp2 = fitlm(pop_ssp2, carbon);

% Predict future carbon emissions
predicted_carbon_ssp2 = predict(model_ssp2, pop_ssp2);

% Calculate 95% prediction intervals for ssp2
[pred_carbon_ssp2, CI_carbon_ssp2] = predict(model_ssp2, pop_ssp2, 'Prediction', 'observation');

% Historical data from 1995 to 2022 (in GtCO2)
historical_years = 1995:2022;
historical_emissions = [29.2, 30.2, 31.9, 30.5, 30.8, 30.9, 30.7, 31.6, 33.5, ...
                        34.0, 34.5, 35.8, 36.1, 36.8, 36.7, 38.5, 39.7, 40.3, ...
                        40.1, 40.7, 41.1, 40.1, 40.6, 41.1, 41.6, 39.3, 41.1, 41.5];

% Projection settings for total carbon footprint
years_projection = 2022:1:2100;
y0_projection = 41.5;  % Start projection from 41.5 GtCO2 in 2022

RCBs = [100, 800];  % Remaining Carbon Budgets in GtCO2 for 83% and 800 Gt
probabilities = [83];  % Probability for 83%

% Initialize figure
figure;
hold on;

% Define custom colors
color2 = [0.8500, 0.3250, 0.0980]; % Red
color3 = [0.9290, 0.6940, 0.1250]; % Yellow
color4 = [0.4940, 0.1840, 0.5560]; % Purple
color5 = [0.3010, 0.7450, 0.9330]; % Cyan for RCB 800 GtCO2


% Plot prediction intervals for the construction carbon footprint
fill([year; flipud(year)], [CI_carbon_ssp2(:,1); flipud(CI_carbon_ssp2(:,2))], color4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot predicted values for the construction carbon footprint
plot(year, predicted_carbon_ssp2, 'Color', color4, 'LineWidth', 1.5);

% Plot actual carbon values for the construction carbon footprint
plot(year, carbon,'LineWidth', 2,'Color', color2);

% Plot the historical data for total carbon footprint
plot(historical_years, historical_emissions, 'k', 'LineWidth', 2);

% Loop through each RCB value (100 and 800 GtCO2)
for i = 1:length(RCBs)
    % Define the initial guess for k and the tolerance for total carbon footprint
    k_initial_guess = 0.1;
    
    % Use fzero to find the correct k that results in the exact cumulative emissions (RCB)
    k = fzero(@(k) trapz(years_projection, y0_projection * exp(-k * (years_projection - 2022))) - RCBs(i), k_initial_guess);
    
    % Compute the path using the calculated k for total carbon footprint
    y_projection = y0_projection * exp(-k * (years_projection - 2022));
    
    % Choose color based on RCB value
    if RCBs(i) == 100
        plot(years_projection, y_projection, 'Color', color3, 'LineWidth', 2);  % Yellow for RCB 100
        text(2050, y_projection(find(years_projection==2050)), ...
            ['1.5 degrees'], ...
            'FontSize', 12, 'Color', color3);
    elseif RCBs(i) == 800
        plot(years_projection, y_projection, 'Color', color5, 'LineWidth', 2);  % Cyan for RCB 800
        text(2050, y_projection(find(years_projection==2050)), ...
            ['2 degrees'], ...
            'FontSize', 12, 'Color', color5);
    end
end

% Add labels and title
%title('Predicted Carbon Footprint and Carbon Budget Pathway for ssp2 (83%)');
xlabel('Year');
ylabel('Carbon Footprint / Total CO_2 (GtCO_2)');

% Create a separate legend box
legend({'Confidence Interval', 'Projected Construction Carbon', 'Historical Construction Carbon', 'Global Historical Carbon Emissions', ...
    'Carbon budget 1.5°C (83%)', 'Carbon budget 2°C (83%)'}, ...
    'Location', 'northeastoutside', 'Box', 'on');

% Customize plot
set(gcf, 'Color', 'w');
grid off;
xlim([1995 2050]);
ylim([0 45]);

a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)


hold off;

set(gcf, 'Position',  [387,275,769,553])


% Save the plot
saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/combined_carbon_footprint_and_budget_bothdegrees_ssp2.svg');

%% carbon budget bar graph
% Define data for the 1.5°C scenario
co2_1_5C = [500, 300, 250, 150, 100];  % CO2 emissions in billion tonnes
likelihood_1_5C = [17, 33, 50, 67, 83];  % Probabilities in percentages
global_emissions_2022 = 41;  % Global emissions in 2022

% Define data for the 2°C scenario
co2_2C = [2000, 1450, 1150, 950, 800];  % CO2 emissions in billion tonnes
likelihood_2C = [17, 33, 50, 67, 83];  % Probabilities in percentages

% Create a figure
figure;
set(gcf, 'Position', [100, 100, 1200, 600]);  % Set the figure size

% Create the subplot for 1.5°C scenario
subplot(1, 2, 1);
bar(likelihood_1_5C, co2_1_5C, 'FaceColor', [0.4, 0.6, 0.8]);  % Light blue bars
ylim([0, 550]);  % Set the y-axis limits
xlim([10, 90]);  % Set the x-axis limits for spacing
set(gca, 'XTick', [likelihood_1_5C, 95], 'XTickLabel', {'17%', '33%', '50%', '67%', '83%'});
ylabel('CO_2 emissions (billion tonnes)');
title('Keeping global average temperature rise below 1.5°C');

% Create the subplot for 2°C scenario
subplot(1, 2, 2);
bar(likelihood_2C, co2_2C, 'FaceColor', [0.4, 0.6, 0.6]);  % Greenish bars
ylim([0, 2100]);  % Set the y-axis limits
xlim([10, 90]);  % Set the x-axis limits for spacing
set(gca, 'XTick', [likelihood_2C, 95], 'XTickLabel', {'17%', '33%', '50%', '67%', '83%'});
title('Keeping global average temperature rise below 2°C');
ylabel('CO_2 emissions (billion tonnes)');

% Add a common title
sgtitle('CO_2 Emissions and Probabilities for Limiting Global Warming');

% Customize appearance
set(gcf, 'Color', 'w');


saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/carbonbudget_bar.svg');
%% construction bar+carbon budget line
construction_carbon_2030=sum(predicted_carbon_ssp2(29:36,:))% from 2023-2030
construction_carbon_2040=sum(predicted_carbon_ssp2(29:46,:))% from 2023-2040
construction_carbon_2050=sum(predicted_carbon_ssp2(29:56,:))% from 2023-2030

bar([construction_carbon_2030,construction_carbon_2040,construction_carbon_2050])

% Add horizontal lines at the specified y-values
yline(100, '--', 'LineWidth', 1.5);   % Line at y = 100
yline(250, '--', 'LineWidth', 1.5);   % Line at y = 250
yline(300, '--', 'LineWidth', 1.5);   % Line at y = 350

a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)

saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/carbonbudget_construction_bar.svg');

% predicted_carbon_ssp3cum=cumsum(predicted_carbon_ssp3(29:56,:))
% plot(predicted_carbon_ssp3cum)
%% Final figure (cumulative+annual)
%% Final figure (cumulative+annual)
%Figure 1
subplot(1,3,1)

construction_carbon_2030=sum(predicted_carbon_ssp2(29:36,:))% from 2023-2030
construction_carbon_2040=sum(predicted_carbon_ssp2(29:46,:))% from 2023-2040
construction_carbon_2050=sum(predicted_carbon_ssp2(29:56,:))% from 2023-2030

bar([construction_carbon_2030,construction_carbon_2040,construction_carbon_2050])
ylabel('Cumulative CO_2 from Construction (Gt)');
xticklabels({'2030','2040','2050'})
xlabel('Year')

% Add horizontal lines at the specified y-values
yline(100, '--', 'LineWidth', 1.5);   % Line at y = 100
yline(250, '--', 'LineWidth', 1.5);   % Line at y = 250
yline(300, '--', 'LineWidth', 1.5);   % Line at y = 350

a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',16)

subplot(1,3,[2,3])
% Figure 2
% Load the Excel file
filename = '/Users/lichaohui/Desktop/calculation/construction/r_forecast/world/regression.xlsx';

% Read the data from the sheet named 'merge'
opts = detectImportOptions(filename, 'Sheet', 'merge');
regression_data = readtable(filename, opts);

% Extract variables for linear regression
year = regression_data.year;
carbon = regression_data.carbon / 1e12; % Convert to GtCO2
pop = regression_data.pop;
pop_ssp2 = regression_data.pop_ssp2;

% Perform linear regression
model_ssp2 = fitlm(pop_ssp2, carbon);

% Predict future carbon emissions
predicted_carbon_ssp2 = predict(model_ssp2, pop_ssp2);

% Calculate 95% prediction intervals for ssp2
[pred_carbon_ssp2, CI_carbon_ssp2] = predict(model_ssp2, pop_ssp2, 'Prediction', 'observation');

% Historical data from 1995 to 2022 (in GtCO2)
historical_years = 1995:2022;
historical_emissions = [29.2, 30.2, 31.9, 30.5, 30.8, 30.9, 30.7, 31.6, 33.5, ...
                        34.0, 34.5, 35.8, 36.1, 36.8, 36.7, 38.5, 39.7, 40.3, ...
                        40.1, 40.7, 41.1, 40.1, 40.6, 41.1, 41.6, 39.3, 41.1, 41.5];

% Projection settings for total carbon footprint
years_projection = 2022:1:2100;
y0_projection = 41.5;  % Start projection from 41.5 GtCO2 in 2022

RCBs = [100, 800];  % Remaining Carbon Budgets in GtCO2 for 83% and 800 Gt
probabilities = [83];  % Probability for 83%

% Initialize figure

% Define custom colors
color2 = [0.8500, 0.3250, 0.0980]; % Red
color3 = [0.9290, 0.6940, 0.1250]; % Yellow
color4 = [0.4940, 0.1840, 0.5560]; % Purple
color5 = [0.3010, 0.7450, 0.9330]; % Cyan for RCB 800 GtCO2


% Plot confidence intervals first, so lines will be on top
fill([year; flipud(year)], [CI_carbon_ssp2(:,1); flipud(CI_carbon_ssp2(:,2))], color4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot predicted values
hold on;
plot(year, predicted_carbon_ssp2, 'Color', color4, 'LineWidth', 1.5);

% Plot actual carbon values for construction
plot(year, carbon,'LineWidth', 2,'Color', color2);

% Plot historical data for total carbon footprint
plot(historical_years, historical_emissions, 'k', 'LineWidth', 2);

% Loop through RCB values and plot them with adjusted text and color
for i = 1:length(RCBs)
    k = fzero(@(k) trapz(years_projection, y0_projection * exp(-k * (years_projection - 2022))) - RCBs(i), k_initial_guess);
    y_projection = y0_projection * exp(-k * (years_projection - 2022));
    
    if RCBs(i) == 100
        plot(years_projection, y_projection, 'Color', color3, 'LineWidth', 2);
        text(2050, y_projection(find(years_projection==2050)), '1.5 degrees', 'FontSize', 12, 'Color', color3);
    elseif RCBs(i) == 800
        plot(years_projection, y_projection, 'Color', color5, 'LineWidth', 2);
        text(2050, y_projection(find(years_projection==2050)), '2 degrees', 'FontSize', 12, 'Color', color5);
    end
end

% Customize the plot appearance
xlabel('Year');
ylabel('Construction CO_2 Footprint and Carbon Budget Trajectory (Gt)');
xlim([1995 2050]);
ylim([0 45]);

legend({'Confidence Interval', 'Projected Construction Carbon', 'Historical Construction Carbon', 'Global Historical Carbon Emissions', ...
        'Carbon budget 1.5°C (83%)', 'Carbon budget 2°C (83%)'}, 'Location', 'northeast', 'Box', 'on');

set(gca, 'FontSize', 16);
hold off;


set(gcf, 'Position',  [121,318,1349,564])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/Final_figure1.svg');

%%

% Given CO2_Wfootprint_constructionTimeseries data
CO2_Wfootprint_constructionTimeseries = [0.5, 0.6, 0.65, 0.7, 0.75, 0.78, 0.8, 0.83, 0.85, 0.88, 0.9, 0.93, 0.96, ...
    1.0, 1.05, 1.08, 1.12, 1.15, 1.2, 1.25, 1.3, 1.35, 1.38, 1.42, 1.45, 1.5, 1.55, 1.6];

% Corresponding years
years = (1995:2022)';

% Convert to time series object
data = CO2_Wfootprint_constructionTimeseries';
ts = array2timetable(data, 'RowTimes', years, 'VariableNames', {'CO2_Wfootprint'});

% Check for stationarity (ADF Test)
[h,pValue,stat,cValue,reg] = adftest(ts.CO2_Wfootprint);

% Fit ARIMA model (1,1,1) since the data is non-stationary
model = arima(1,1,1);
fit = estimate(model, ts.CO2_Wfootprint);

% Forecast future values until 2050
numYearsAhead = 2050 - 2022;
[forecastedValues, forecastMSE] = forecast_try(fit, numYearsAhead, 'Y0', ts.CO2_Wfootprint);

% Create future years for plotting
futureYears = (2023:2050)';

% Plot the historical data and the forecasted values
figure;
plot(years, ts.CO2_Wfootprint, 'b', 'LineWidth', 2); hold on;
plot(futureYears, forecastedValues, 'r--', 'LineWidth', 2);
xlabel('Year');
ylabel('CO2 Wfootprint');
title('CO2 Wfootprint Construction Timeseries with ARIMA Forecast');
legend('Historical Data', 'Forecasted Data');
grid on;

%%
% Your historical data (CO2_Wfootprint_constructionTimeseries)
years = 1995:2022;

% Fit ARIMA model
model = arima(1,1,0); % You can modify the order as needed
fitModel = estimate(model, CO2_Wfootprint_constructionTimeseries');

% Forecast future values until 2050
future_years = 2023:2050;
num_future_years = length(future_years);
[forecasted_values, ~] = forecast(fitModel, num_future_years, 'Y0', CO2_Wfootprint_constructionTimeseries');

% Plot historical and projected values
figure;
plot(years, CO2_Wfootprint_constructionTimeseries, '-o', 'DisplayName', 'Historical Data');
hold on;
plot(future_years, forecasted_values, '-x', 'DisplayName', 'Projected Data');
xlabel('Year');
ylabel('CO2 Footprint');
title('ARIMA Projected Construction CO2 Footprint');
legend;
grid on;
hold off;

saveas(gcf,'construction/figures/arima.svg')

