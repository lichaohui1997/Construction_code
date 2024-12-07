% By Chaohui Li
clear;
clc;
Path = '/Users/lichaohui/Desktop/calculation';
%% Set working directory
cd(Path)
%% Useful functions
addpath('command/paneldata');
addpath('command/paneldata/util');
addpath('command/paneldata/tests');
addpath('command/paneldata/stafun');
addpath('command/paneldata/estimation');
%% Load the Y variable
load('construction/regression_Y_variable.mat')

Bluewater_Wfootprint_regionsTimeseries=Bluewater_Wfootprint_regionsTimeseries(:,1:27);
CO2_Wfootprint_regionsTimeseries=CO2_Wfootprint_regionsTimeseries(:,1:27);
Energy_Wfootprint_regionsTimeseries=Energy_Wfootprint_regionsTimeseries(:,1:27);

Energy_Wfootprint_constructionTimeseries=Energy_Wfootprint_constructionTimeseries(:,1:27);
CO2_Wfootprint_constructionTimeseries=CO2_Wfootprint_constructionTimeseries(:,1:27);
Bluewater_Wfootprint_constructionTimeseries=Bluewater_Wfootprint_constructionTimeseries(:,1:27);

gEnergy=Energy_Wfootprint_constructionTimeseries;
gCO2=CO2_Wfootprint_constructionTimeseries;
gBluewater=Bluewater_Wfootprint_constructionTimeseries
%%
CO2_1=CO2_Wfootprint_regionsTimeseries;
CO2_1(45:49,:)=[];
CO2_1(41,:)=[];
CO2_tot=reshape(CO2_1',numel(CO2_1),1);
CO2_tot1=reshape(CO2_tot,numcountries,numyear)
CO2_per=CO2_tot./population
CO2_per1=reshape(CO2_per,numcountries,numyear)
%
Bluewater_1=Bluewater_Wfootprint_regionsTimeseries;
Bluewater_1(45:49,:)=[];
Bluewater_1(41,:)=[];
Bluewater_tot=reshape(Bluewater_1',numel(Bluewater_1),1);
Bluewater_per=Bluewater_tot./population
Bluewater_per1=reshape(Bluewater_per,numcountries,numyear)
%
Energy_1=Energy_Wfootprint_regionsTimeseries;
Energy_1(45:49,:)=[];
Energy_1(41,:)=[];
Energy_tot=reshape(Energy_1',numel(Energy_1),1);
Energy_per=Energy_tot./population
Energy_per1=reshape(Energy_per,numcountries,numyear)

%% Start the regression process
%% CO2

Country = reshape(repmat(countries, numyear,1),[],1);
Years_reg = reshape(repmat(1995:2021, numcountries,1)',[],1);

Country = Country(2:end);
Years_reg = Years_reg(2:end);

y1=diff(log(abs(CO2_tot) + 1e-6));
y2=diff(log(abs(Energy_tot) + 1e-6));
y3=diff(log(abs(Bluewater_tot) + 1e-6));


x1=diff(log(abs(gdp) + 1e-6));
x2=diff(log(abs(industryva) + 1e-6));
x3=diff(log(urbanpopulation(:)));
x4=diff(log(abs(population) + 1e-6));
x5=log(abs(gdp) + 1e-6).^2
x5(1,:)=[]


Years_reg=num2cell(Years_reg) % this is important, to make the year into string cell, so that fixed effect model can be carried out using year dummy
Years_reg = cellfun(@num2str, Years_reg, 'UniformOutput', false);

data = table(Country, Years_reg, y1, y2, y3, x1, x2, x3, x4, x5,...
    'VariableNames', {'Country', 'Year', 'CO2', 'Energy','Bluewater','GDP','Industry','Urban' ,'Population','GDPsq'});

% panel regression models
model1 = fitlm(data, 'CO2 ~ GDP + Industry + Urban + Population + GDPsq');
model2 = fitlm(data, 'CO2 ~ GDP');
model3 = fitlm(data, 'CO2 ~ Population');
model4 = fitlm(data, 'CO2 ~ Urban');
model5 = fitlm(data, 'CO2 ~ Population + Industry');
model6 = fitlm(data, 'CO2 ~ GDP + Population');

% fixed effect models
model7 = fitlm(data, 'CO2 ~ Population + Country + Year');


% Display regression coefficients and statistics
for i = 1:7
  coeffs_table{i} = eval(sprintf('model%d.Coefficients', i)); % Store the coefficients in the cell array
  pvals = coeffs_table{i}.pValue;
  sig_star = cell(size(pvals));
  sig_star(pvals <= 0.1) = {'*'};
  sig_star(pvals <= 0.05) = {'**'};
  sig_star(pvals <= 0.01) = {'***'};
  coeffs_table{i}.sig_star = sig_star;   
    
    coeffs_table{i} = [coeffs_table{i}];
end

for i = 1:length(coeffs_table)
    writetable(coeffs_table{i}, '/Users/lichaohui/Desktop/calculation/construction/coeffs_table.xlsx', 'Sheet', i, 'WriteRowNames', true);
end

%%
%% KPSS test for stationarity. Result: Stationary
[h_kpss, p_kpss] = kpsstest(y1);
if h_kpss == 1
    disp('The data is not stationary according to the KPSS test');
else
    disp('The data is stationary according to the KPSS test');
end
%% Loop KPSS test for stationarity for the 9 variables
% Load data and define variables
varNames = {'x1', 'x2', 'x3', 'x4','y1', 'y2', 'y3'};
nVars = length(varNames);

% Create table to store results
kpssTable = table('Size',[nVars 3],'VariableTypes',{'string','double','double'},'VariableNames',{'Variable','TestStatistic','pValue'});

% Perform KPSS tests for each variable and store results in table
for i = 1:nVars
    % Perform KPSS test
    [h,pValue,testStat] = kpsstest(eval(varNames{i}));
    
    % Store results in table
    kpssTable(i,:) = {varNames{i}, testStat, pValue};
end

% Display table of KPSS test results
disp(kpssTable)

writetable(kpssTable,'/Users/lichaohui/Desktop/calculation/construction/manuscript/revised/DataS3.xlsx','Sheet','kpss_test')

%% Test for stationarity with PP test. Result: no unit root
[h, pValue, stat, cValue, reg] = pptest(y1, 'model', 'ar', 'lags', 4);
if h == 1
    disp('The null hypothesis of a unit root is rejected at the 5% significance level.')
else
    disp('The null hypothesis of a unit root cannot be rejected at the 5% significance level.')
end
%% Loop PP test for 9 variables
% Create table to store results
ppTable = table('Size',[nVars 4],'VariableTypes',{'string','double','double','double'},'VariableNames',{'Variable','pValue','Statistics','cValue'});

% Perform KPSS tests for each variable and store results in table
for i = 1:nVars
    % Perform KPSS test
    [h, pValue, stat, cValue, reg] = pptest(y1, 'model', 'ar', 'lags', 4);

    % Store results in table
    ppTable(i,:) = {varNames{i}, pValue,stat,cValue};
end

% Display table of KPSS test results
disp(ppTable)

writetable(ppTable,'/Users/lichaohui/Desktop/calculation/construction/manuscript/revised/DataS3.xlsx','Sheet','pp_test')

%% % Perform ADF test for stationarity. Result: Not stationary
adf_test = adftest(y1);
disp('ADF test:');
disp(adf_test);
%
adf_test = adftest(x1);
disp('ADF test:');
disp(adf_test);
%% Correlation test.you can continue with panel analysis if the Pearson correlation coefficient is statistically significant
% Calculate the Pearson correlation coefficient
% r = corrcoef(x1, y);
% 
% % Display the result
% disp(['The Pearson correlation coefficient between x and y is ', num2str(r(1, 2))]);
% 
% % x and y are your variables
% % alpha is the significance level
% 
% [r,p] = corr(x1,y); % calculate the correlation coefficient and p-value
% n = length(x1); % sample size
% t = r * sqrt(n - 2) / sqrt(1 - r^2); % calculate t-statistic
% 
% df = n - 2; % degrees of freedom
% alpha = 0.05
% tcrit = tinv(1 - alpha/2, df); % critical t-value
% 
% if abs(t) > tcrit
%     fprintf('The correlation coefficient is statistically significant.\n');
% else
%     fprintf('The correlation coefficient is not statistically significant.\n');
% end
%% Ljung-Box test for autocorrelation. Result: No autocorrelation
max_lag = 10; % Set the maximum number of lags
alpha = 0.05; % Set the significance level

% Compute the autocorrelation function (ACF) of the residuals
residuals = model7.Residuals.Raw;
[acf, lags, bounds] = autocorr(residuals, max_lag);

% Compute the Ljung-Box test statistic and p-value
[LB_stat, LB_pvalue] = lbqtest(residuals, 'Lags', 1:max_lag, 'Alpha', alpha);

% Display the results
disp('Ljung-Box Test for Autocorrelation:');
disp(['Maximum Lag: ', num2str(max_lag)]);
disp(['Significance Level: ', num2str(alpha)]);
disp(' ');
disp('Autocorrelation Function (ACF):');
%disp([lags', acf, bounds]);

disp(' ');
disp(['Ljung-Box Test Statistic: ', num2str(LB_stat)]);
disp(['P-value: ', num2str(LB_pvalue)]);
if LB_pvalue < alpha
    disp('Autocorrelation detected (Ljung-Box test)');
else
    disp('No autocorrelation detected (Ljung-Box test)');
end
%% Ljung-Box test for autocorrelation for 10 models
% Initialize cell array to store test results
test_results = cell(1,6);

% Loop through the 10 models
for i = 1:6
    % Get the residuals from the model
    residuals = eval(sprintf('model%d.Residuals.Raw', i));
    
    % Perform the Ljung-Box test with 20 lags
    [h, p] = lbqtest(residuals, 'Lags', 20);
    
    % Store the test results in a table
    test_results{i} = table(h, p, 'VariableNames', {'Ljung-Box_Stat', 'P_Value'});
end

% Combine the test results into one table
test_results_table = vertcat(test_results{:});

writetable(test_results_table,'/Users/lichaohui/Desktop/calculation/construction/manuscript/revised/DataS3.xlsx','Sheet','Ljung_Box_test')

%% Durbin-Watson test for autocorrelation. Result: Autocorrelation detected, because the result number is not close to 2 (range:0-4)
autocorrelation_test = dwtest(model1);
disp('Autocorrelation Test:');
disp(autocorrelation_test);
%%
%% U-shape residual check
fitted = predict(model1);
residuals = model1.Residuals.Raw;
scatter(fitted, residuals);
xlabel('Fitted Values');
ylabel('Residuals');
title('U-Shaped Relationship Check');

%% Display table of residuals
residuals = table(model7.Residuals.Raw, model7.Fitted);
residuals.Properties.VariableNames = {'Residuals', 'Fitted Values'};
disp('Residuals:');
disp(residuals);

% Plot residuals
plotResiduals(model7, 'fitted');
title('Residuals', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Fitted Values', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Residuals', 'FontSize', 12, 'FontWeight', 'bold');

%% Compute correlation coefficients
corr_coeffs = corrcoef(data{:, 3:end});

% Plot heatmap of correlation coefficients
figure;
heatmap(corr_coeffs, 'Colormap', colormap('redbluecmap'), 'ColorbarVisible', 'off', 'FontSize', 10);
title('Correlation Matrix');
xlabel('Variables');
ylabel('Variables');


%% heteroscedasticity check. It will show no heteroscedasticity detected for all models but it's because pValue=NaN

[h, pValue, stat, cv] = archtest(model1.Residuals.Raw, 'Lags', 1:5);

% Display results
if pValue < 0.05
    disp('Heteroscedasticity detected (Goldfeld-Quandt test)');
else
    disp('No heteroscedasticity detected (Goldfeld-Quandt test)');
end
%% heteroscedasticity check. GARCH model
Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
EstMdl = estimate(Mdl,x1);

%% garch for 10 models

%% heteroscedasticity check. GARCH model
% Initialize cell array to store model estimates
model_estimates = cell(1,6);

% Loop through the 10 models
for i = 1:6
    % Get the residuals from the model
    residuals = eval(sprintf('model%d.Residuals.Raw', i));
    
    % Specify the GARCH model
    Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
    
    % Estimate the GARCH model parameters
    EstMdl = estimate(Mdl, residuals);
    
    % Store the estimated models
    model_estimates{i} = EstMdl;
end
%% present garch result model
%% Create table of GARCH model estimates

%% GARCH models and heteroskedasticity checks for 10 models
% Initialize cell array to store model estimates
model_estimates = cell(1,6);

% Loop through the 10 models
for i = 1:6
    % Get the residuals from the model
    residuals = eval(sprintf('model%d.Residuals.Raw', i));
    
    % Define a GARCH(1,1) model
    Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
    
    % Estimate the GARCH model with the residuals
    EstMdl = estimate(Mdl, residuals);
    
    % Store the model estimates in the cell array
    model_estimates{i} = EstMdl;
end

%% Create table of GARCH model estimates and statistics

%% GARCH models and heteroskedasticity checks for 10 models
% Initialize cell array to store model estimates
model_estimates = cell(1,6);

% Loop through the 10 models
for i = 1:6
    % Get the residuals from the model
    residuals = eval(sprintf('model%d.Residuals.Raw', i));
    
    % Define a GARCH(1,1) model
    Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);
    
    % Estimate the GARCH model with the residuals
    EstMdl = estimate(Mdl, residuals);
    
    % Store the model estimates in the cell array
    model_estimates{i} = EstMdl;
end


%% Create table of GARCH model estimates
%% Create table of GARCH model estimates
% Initialize arrays to store estimates
P = zeros(6,1);
Q = zeros(6,1);
Constant = zeros(6,1);
Offset = zeros(6,1);
GARCH = zeros(6,1);
ARCH = zeros(6,1);
UnconditionalVariance = zeros(6,1);

% Loop through the 10 models
for i = 1:6
    % Get the estimates
    P(i) = model_estimates{i}.P;
    Q(i) = model_estimates{i}.Q;
    Constant(i) = model_estimates{i}.Constant;
    Offset(i) = model_estimates{i}.Offset;
    GARCH(i) = model_estimates{i}.GARCH{1}; % extract value from cell
    ARCH(i) = model_estimates{i}.ARCH{1}; % extract value from cell
    UnconditionalVariance(i) = model_estimates{i}.UnconditionalVariance;
end

% Create table of estimates
EstimatesTable = table(P, Q, Constant, Offset, GARCH, ARCH, UnconditionalVariance, ...
    'VariableNames', {'P', 'Q', 'Constant', 'Offset', 'GARCH', 'ARCH', 'UnconditionalVariance'}, ...
    'RowNames', {'Model1', 'Model2', 'Model3', 'Model4', 'Model5', 'Model6'});

% Display the table
disp('GARCH Model Estimates:');
disp(EstimatesTable);


writetable(EstimatesTable,'/Users/lichaohui/Desktop/calculation/construction/manuscript/revised/DataS3.xlsx','Sheet','GARCH_test')

%% Trying to create presentatble graphs
vec = repmat(1:17, numcountries, 1)'
vec=vec(:);
%scatter(urbanpopulation,CO2_tot,[],vec)

