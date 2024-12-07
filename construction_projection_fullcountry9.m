%% Population
% Load the Excel data into MATLAB
pop_1 = readcell('/Users/lichaohui/Desktop/calculation/construction/P_Data_Extract_From_Population_estimates_and_projections_regions.xlsx', 'Sheet', 'Sheet1', 'Range', 'A1:BE258');

% Define regions as cell arrays of country codes
India = {'IND'};
China = {'CHN'};

NorthAmerica = {'ABW', 'ATG', 'PYF', 'CAN', 'USA', 'MEX', 'BLZ', 'GTM', 'HND', 'SLV', 'NIC', 'CRI', 'PAN', 'PRI', 'VIR', 'BHS', ...
  'BRB', 'DMA', 'GRD', 'KNA', 'LCA', 'VCT', 'TTO', 'VGB', 'CUW', 'SXM', 'CYM', 'BMU'};

EU = {'CHE', 'GBR', 'ISL', 'MCO', 'NOR', 'SMR', 'AUT', 'BEL', 'BGR', 'HRV', 'CYP', 'CZE', 'DNK', 'EST', 'FIN', 'FRA', 'DEU', ...
  'GRC', 'HUN', 'IRL', 'ITA', 'LVA', 'LIE', 'LTU', 'LUX', 'MLT', 'NLD', 'POL', 'PRT', 'ROU', 'SVK', 'SVN', 'ESP', 'SWE'};

Africa = {'BWA', 'SOM', 'TUN', 'AGO', 'DZA', 'BEN', 'BFA', 'BDI', 'CPV', 'CMR', 'CAF', 'TCD', 'COM', 'COD', 'COG', 'CIV', ...
  'DJI', 'EGY', 'GNQ', 'ERI', 'ETH', 'GAB', 'GMB', 'GHA', 'GIN', 'GNB', 'KEN', 'LSO', 'LBR', 'LBY', 'MDG', 'MWI', ...
  'MLI', 'MRT', 'MUS', 'MAR', 'MOZ', 'NAM', 'NER', 'NGA', 'RWA', 'STP', 'SEN', 'SYC', 'SLE', 'ZAF', 'SSD', 'SDN', ...
  'SWZ', 'TGO', 'TZA', 'UGA', 'ZMB', 'ZWE'};

OtherAsiaOceania = {'BGD', 'MDV', 'MYS', 'NCL', 'NPL', 'PAK', 'THA', 'TUR', 'AFG', 'ASM', 'AUS', 'BRN', 'BTN', 'KHM', ...
  'FJI', 'FSM', 'HKG', 'IDN', 'JPN', 'KIR', 'LAO', 'MAC', 'MHL', 'MNG', 'MMR', 'NRU', 'NZL', 'PLW', ...
  'PNG', 'PHL', 'KOR', 'PRK', 'SGP', 'SLB', 'WSM', 'TLS', 'TON', 'TUV', 'VUT', 'VNM'};

OtherEurope = {'ALB', 'ARM', 'AZE', 'BLR', 'BIH', 'GEO', 'KAZ', 'XKX', 'KGZ', 'MDA', 'MNE', 'MKD', 'RUS', 'SRB', 'TJK', ...
  'TKM', 'UKR', 'UZB', 'AND', 'GIB', 'IMN'};

LatinAmerica = {'ARG', 'BOL', 'BRA', 'CHL', 'COL', 'CUB', 'DOM', 'ECU', 'GUY', 'HTI', 'JAM', 'PER', 'PRY', 'URY', 'VEN', 'SUR'};

MiddleEast = {'BHR', 'IRN', 'IRQ', 'ISR', 'JOR', 'KWT', 'LBN', 'OMN', 'QAT', 'SAU', 'SYR', 'ARE', 'PSE', 'YEM'};

OtherCountries = {'ARB', 'CSS', 'CHI', 'FCS', 'FRO', 'GUM', 'GRL', 'MAF', 'MNP', 'OSS', 'PSS', 'TCA', 'WLD', 'IBD', ...
  'IBT', 'IDB', 'IDX', 'IDA', 'INX', 'MIC', 'LIC', 'LMC', 'LMY', 'HIC', 'UMC', 'IBRD only', 'OED', 'OSS', ...
  'HPC', 'EAR', 'LTE', 'PRE', 'PST', 'NAC', 'HPC', 'HIC', 'HIPC', 'SSA', 'LMY', 'UMC', 'HIC'};

% Combine all regions into a cell array
Regions = {India, China, NorthAmerica, EU, Africa, OtherAsiaOceania, OtherEurope, LatinAmerica, MiddleEast, OtherCountries};
RegionNames = {'India', 'China', 'NorthAmerica', 'EU', 'Africa', 'OtherAsiaOceania', 'OtherEurope', 'LatinAmerica', 'MiddleEast', 'OtherCountries'};

% Extract the first column containing country codes from the Excel data
CountryOrder = pop_1(:, 1); % Assuming the first column contains the country codes

% Initialize an empty matrix to store the summed populations
numColumns = size(pop_1, 2) - 1; % Exclude the country code column
pop_2 = zeros(length(Regions), numColumns); % Rows: regions, Columns: years

% Loop through each region and sum the population for each year
for i = 1:length(Regions)
region = Regions{i};
regionData = []; % To store data for countries in the region
% Loop through the region's country codes and extract the relevant rows from pop_1
    for j = 1:length(region)
        idx = find(strcmp(CountryOrder, region{j}));
        if ~isempty(idx)
            regionData = [regionData; cell2mat(pop_1(idx, 2:end))]; % Collect population data for that country
        end
    end
    % Sum the population data across countries for each year
    if ~isempty(regionData)
        pop_2(i, :) = sum(regionData, 1);
    end
end

% Display the result (pop_2)
disp(array2table(pop_2, 'RowNames', RegionNames));


pop_2([10], :) = []

pop_proj_region_sum=pop_2

%%
load('regression_Y_variable.mat')
writematrix(CO2_Wfootprint_regionsTimeseries,'/Users/lichaohui/Desktop/calculation/construction/regression_Y_variable.xls')

var_y=CO2_Wfootprint_regionsTimeseries


Region9 = {"India", "China", "NorthAmerica", "EU", "Africa", ...
           "OtherAsiaOceania", "OtherEurope", "LatinAmerica", "MiddleEast"};


% Assuming CO2_Wfootprint_regionsTimeseries is a matrix of size 41x28
var_y1 = mat2cell(var_y, ones(49,1), ones(28,1));


India = [35];                 % India
China = [31,41];                 % China
NorthAmerica=[32,29]
EU = [1:28];                  % Austria, Belgium, Bulgaria, Cyprus, Czech Republic, Germany, Denmark, Estonia, Spain, Finland, France, Greece, Croatia, Hungary, Ireland, Italy, Lithuania, Luxembourg, Latvia, Malta, Netherlands, Poland, Portugal, Romania, Sweden, Slovenia, Slovakia, United Kingdom
Africa = [44, 48];            % South Africa, RoW Africa
OtherAsiaOceania = [30,33,37,38,43,45];         % Australia, Indonesia
OtherEurope = [39, 40, 42, 47];  % Switzerland, Turkey, Norway, RoW Europe
LatinAmerica = [34,36,46];      % Mexico, RoW America
MiddleEast = [49];            % RoW Middle East

RegionOrder = [India, China, NorthAmerica, EU, Africa, OtherAsiaOceania, OtherEurope, LatinAmerica, MiddleEast]';

aggregationRegion=[1,2,2,28,2,6,4,3,1];
NNNR = 9; %New number of regions


m=1;
var_y2 = cell(size(var_y1));
for k=1:49
for i=1:49
if RegionOrder(k,1)==i
var_y2(m,:)=var_y1(i,:);
m=m+1;
end
end
end

var_y3=mat2cell(var_y2,aggregationRegion,ones(28,1))
var_y4 = cellfun(@(x) sum(cell2mat(cellfun(@(y) sum(y(:), 'omitnan') * isnumeric(y), x, 'UniformOutput', false))), var_y3, 'UniformOutput', false);

var_y5=cell2mat(var_y4)
var_y6=[var_y5,NaN(9,28)]
var_y7=reshape(var_y6,[],1)



%% Yearid and Regionid
year=reshape(repmat(1995:1:2050,9,1),[],1)
yearid=categorical(reshape(repmat(1:1:56,9,1),[],1));

regionid_1=repmat(Region9,1,56)
regionid_2=reshape(regionid_1,[],1)

regionid_3=1:1:9
regionid_4=repmat(regionid_3,1,56)
regionid_5=reshape(regionid_4,[],1)
% establish table
regionid=categorical(regionid_5)
regionname=regionid_2

id=categorical([1:504]');

x=reshape(pop_proj_region_sum,[],1)
% x=diff(log(abs(x)+ 1e-6));
% y=diff(log(abs(var_y7)+ 1e-6));

y=var_y7;

% id=id(2:end)
% yearid=yearid(2:end)
% regionid=regionid(2:end)

RegTable=table(id,yearid,regionid,x,y)
writetable(RegTable,'/Users/lichaohui/Desktop/calculation/construction/RegTable.csv')

%% start regression

% Load the data
RegTable = readtable('/Users/lichaohui/Desktop/calculation/construction/RegTable.csv');

% Convert yearid and regionid to categorical types (equivalent to factors in R)
RegTable.yearid = categorical(RegTable.yearid);
RegTable.regionid = categorical(RegTable.regionid);

% Get the list of unique regions
unique_regions = unique(RegTable.regionid);

% Loop through each region and perform OLS regression
for i = 1:length(unique_regions)
    
    % Filter data for the current region
    region_data = RegTable(RegTable.regionid == unique_regions(i), :);
    
    % Perform linear interpolation for missing population data
    % In MATLAB, we use `fillmissing` with 'linear' for interpolation
    region_data.x = fillmissing(region_data.x, 'linear');
    
    % Perform OLS regression: y ~ x (carbon footprint ~ population)
    % Linear model for OLS regression
    mdl = fitlm(region_data.x, region_data.y);
    
    % Find future data (years 29-56)
    future_data = region_data(str2double(string(region_data.yearid)) > 28, :);
    
    % Predict future carbon emissions
    future_data.predicted_carbon = predict(mdl, future_data.x);
    
    % 95% prediction intervals for future data
    [predicted_carbon, CI] = predict(mdl, future_data.x, 'Alpha', 0.05);
    future_data.lower_bound = CI(:, 1);  % Lower bound of 95% CI
    future_data.upper_bound = CI(:, 2);  % Upper bound of 95% CI
    
    % Predict historical values for y based on the fitted model
    region_data.predicted_carbon = predict(mdl, region_data.x);
    
    % Create the plot for the current region
    figure;
    hold on;
    
    % Plot the 95% confidence interval for the future data
    fill([future_data.yearid; flipud(future_data.yearid)], ...
         [future_data.lower_bound; flipud(future_data.upper_bound)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
     
    % Plot the predicted carbon footprint (historical + future)
    plot(region_data.yearid, region_data.predicted_carbon, 'r-', 'LineWidth', 2);
    
    % Plot the actual historical carbon footprint
    plot(region_data.yearid, region_data.y, 'k--', 'LineWidth', 2);
    
    % Set title and labels
    title(['Carbon Footprint Prediction for Region ', char(unique_regions(i))]);
    xlabel('Year');
    ylabel('Carbon Footprint');
    
    hold off;
    
    % Optionally save the plot
    % saveas(gcf, ['carbon_footprint_region_', char(unique_regions(i)), '.pdf']);
end
%% put them all in one plot
% Load the data
RegTable = readtable('/Users/lichaohui/Desktop/calculation/construction/RegTable.csv');

% Convert yearid and regionid to categorical types (equivalent to factors in R)
RegTable.yearid = categorical(RegTable.yearid);
RegTable.regionid = categorical(RegTable.regionid);

% Get the list of unique regions
unique_regions = unique(RegTable.regionid);

% Define the region names in the correct order
Region9 = {"India","China",  "NorthAmerica", "EU", "Africa",  ...
           "OtherAsiaOceania", "OtherEurope", "LatinAmerica", "MiddleEast"}';
               
               % Create a tiled layout for the plots (assuming 14 regions)
               num_regions = length(unique_regions);
               tiledlayout(4, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % Adjust based on the number of regions and plot size
               
               % Loop through each region and perform OLS regression
               for i = 1:num_regions
               
               % Filter data for the current region
               region_data = RegTable(RegTable.regionid == unique_regions(i), :);
               
               % Perform linear interpolation for missing population data
               region_data.x = fillmissing(region_data.x, 'linear');
               
               % Perform OLS regression: y ~ x (carbon footprint ~ population)
               mdl = fitlm(region_data.x, region_data.y);
               
               % Find future data (years 29-56)
               future_data = region_data(str2double(string(region_data.yearid)) > 28, :);
               
               % Predict future carbon emissions
               future_data.predicted_carbon = predict(mdl, future_data.x);
               
               % 95% prediction intervals for future data
               [predicted_carbon, CI] = predict(mdl, future_data.x, 'Alpha', 0.05);
               future_data.lower_bound = CI(:, 1);  % Lower bound of 95% CI
               future_data.upper_bound = CI(:, 2);  % Upper bound of 95% CI
               
               % Predict historical values for y based on the fitted model
               region_data.predicted_carbon = predict(mdl, region_data.x);
               
               % Create the plot for the current region
               nexttile;  % Move to the next subplot
               
               hold on;
               
               % Plot the 95% confidence interval for the future data
               fill([future_data.yearid; flipud(future_data.yearid)], ...
                    [future_data.lower_bound; flipud(future_data.upper_bound)], ...
                    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
               
               % Plot the predicted carbon footprint (historical + future)
               plot(region_data.yearid, region_data.predicted_carbon, 'r-', 'LineWidth', 2);
               
               % Plot the actual historical carbon footprint
               plot(region_data.yearid, region_data.y, 'k--', 'LineWidth', 2);
               
               % Set title and labels
               title(Region9{i}, 'FontSize', 8);  % Use the corresponding region name from Region9
               xlabel('Year', 'FontSize', 6);
               ylabel('Carbon Footprint', 'FontSize', 6);
               
               hold off;
               end
               
               % Optionally, you can save the entire figure
               % saveas(gcf, 'carbon_footprint_all_regions.pdf');
%% improved
% Load the data
RegTable = readtable('/Users/lichaohui/Desktop/calculation/construction/RegTable.csv');
               
% Convert yearid and regionid to categorical types (equivalent to factors in R)
RegTable.yearid = categorical(RegTable.yearid);
RegTable.regionid = categorical(RegTable.regionid);
RegTable.y=RegTable.y/1e12               
% Get the list of unique regions
unique_regions = unique(RegTable.regionid);
               
% Define the region names in the correct order
               
Region9_title  = {"India","China",  "North America", "EU", "Africa",  ...
                 "Other Asia&Oceania", "Other Europe", "Latin America", "Middle East"}';


% Convert yearid to actual years (1995-2050)
% Assume that yearid starts at 1 for 1995 and ends at 56 for 2050
years = linspace(1995, 2050, 56);

% Create a tiled layout for the plots (assuming 14 regions)
num_regions = length(unique_regions);
tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact'); % Adjust based on the number of regions and plot size

% Loop through each region and perform OLS regression
for i = 1:num_regions
    
    % Filter data for the current region
    region_data = RegTable(RegTable.regionid == unique_regions(i), :);
    
    % Perform linear interpolation for missing population data
    region_data.x = fillmissing(region_data.x, 'linear');
    
    % Perform OLS regression: y ~ x (carbon footprint ~ population)
    mdl = fitlm(region_data.x, region_data.y);
    
    % Find future data (years 29-56)
    future_data = region_data(str2double(string(region_data.yearid)) > 28, :);
    
    % Predict future carbon emissions
    future_data.predicted_carbon = predict(mdl, future_data.x);
    
    % 95% prediction intervals for future data
    [predicted_carbon, CI] = predict(mdl, future_data.x, 'Alpha', 0.05);
    future_data.lower_bound = CI(:, 1);  % Lower bound of 95% CI
    future_data.upper_bound = CI(:, 2);  % Upper bound of 95% CI
    
    % Predict historical values for y based on the fitted model
    region_data.predicted_carbon = predict(mdl, region_data.x);
    
    % Create the plot for the current region
    nexttile;  % Move to the next subplot
    
    hold on;
    
    % Plot the 95% confidence interval for the future data
    fill([years(29:end)'; flipud(years(29:end)')], ...
         [future_data.lower_bound; flipud(future_data.upper_bound)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
     
    % Plot the predicted carbon footprint (historical + future)
    plot(years, region_data.predicted_carbon, 'r-', 'LineWidth', 2);
    
    % Plot the actual historical carbon footprint
    plot(years, region_data.y, 'k--', 'LineWidth', 2);
    
    % Set title and labels
    title(Region9_title{i}, 'FontSize', 14);  % Use the corresponding region name from Region9
    
    % Standardize y-tick format to 10^11 for all subplots
    
    % Set custom ylim for Africa and US
  
    if strcmp(Region9_title{i}, 'EU')
        ylim([0 2]); 
    elseif strcmp(Region9_title{i}, 'North America')
        ylim([0.3 1.5]); 
    elseif strcmp(Region9_title{i}, 'Other Europe')
        ylim([0 0.3]); 
    end

    xlim([1995 2050])
ax = gca
ax.FontSize = 14;
    
hold off;
    
end



set(gcf, 'Position',  [225,230,955,747])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/fullregion_projection.svg');

