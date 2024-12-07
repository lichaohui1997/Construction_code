load('regression_Y_variable.mat')

pop_proj=readcell('P_Data_Extract_From_Population_estimates_and_projections_regions.xlsx')

Canada = {'CAN'};
China = {'CHN'};
India = {'IND'};
SouthKorea = {'KOR'};
Russia = {'RUS'};
Japan = {'JPN'};
US = {'USA'};

EU = {'ALB', 'AND', 'AUT', 'BEL', 'BIH', 'BGR', 'HRV', 'CYP', 'CZE', 'DNK', 'EST', ...
      'FIN', 'FRA', 'DEU', 'GRC', 'HUN', 'ISL', 'IRL', 'ITA', 'LVA', 'LTU', 'LUX', ...
      'MLT', 'MDA', 'MNE', 'NLD', 'MKD', 'NOR', 'POL', 'PRT', 'ROU', 'SRB', 'SVK', ...
      'SVN', 'ESP', 'SWE', 'CHE', 'UKR', 'GBR'};
      
Africa = {'DZA', 'AGO', 'BEN', 'BWA', 'BFA', 'BDI', 'CPV', 'CMR', 'CAF', 'TCD', 'COM', ...
          'COD', 'COG', 'CIV', 'DJI', 'EGY', 'GNQ', 'ERI', 'SWZ', 'ETH', 'GAB', 'GMB', ...
          'GHA', 'GIN', 'GNB', 'KEN', 'LSO', 'LBR', 'LBY', 'MDG', 'MWI', 'MLI', 'MRT', ...
          'MUS', 'MOZ', 'NAM', 'NER', 'NGA', 'RWA', 'STP', 'SEN', 'SYC', 'SLE', 'SOM', ...
          'ZAF', 'SSD', 'SDN', 'TGO', 'TUN', 'UGA', 'ZMB', 'ZWE'};
          
Brazil = {'BRA'};

OtherAsiaOceania = {'AFG', 'BGD', 'KHM', 'IDN', 'IRN', 'IRQ', 'ISR', 'JOR', 'KAZ', 'KWT', 'KGZ', ...
             'LAO', 'LBN', 'MAC', 'MDV', 'MNG', 'MMR', 'NPL', 'PRK', 'OMN', 'PAK', 'PHL', ...
             'QAT', 'SAU', 'SGP', 'KWT', 'SYR', 'TJK', 'THA', 'TLS', 'TKM', 'UZB', 'VNM', ...
             'YEM','AUS', 'FJI', 'FSM', 'KIR', 'MHL', 'NRU', 'NZL', 'PLW', 'PNG', 'SLB', ...
           'TON', 'TUV', 'VUT'};
             
OtherEurope = {'ASM', 'ARM', 'AZE', 'BLR', 'GEO', 'KOS', 'RUS'};

LatinAmerica = {'ATG', 'ARG', 'BRB', 'BLZ', 'BOL', 'CHL', 'COL', 'CRI', 'CUB', 'DMA', ...
                'DOM', 'ECU', 'SLV', 'GTM', 'GUY', 'HTI', 'HND', 'JAM', 'MEX', 'NIC', ...
                'PAN', 'PRY', 'PER', 'PRI', 'KNA', 'LCA', 'MAF', 'VCT', 'SUR', 'TTO', ...
                'URY', 'VEN'};
                
           
MiddleEast = {'BHR', 'IRN', 'IRQ', 'ISR', 'JOR', 'KWT', 'LBN', 'OMN', 'QAT', 'SAU', ...
              'SYR', 'ARE', 'YEM'};

OtherCountries = {'ABW', 'CSS', 'CYM', 'CEB', 'CHI', 'CUW', 'EAR', 'EAS', 'EAP', ...
                  'TEA', 'EMU', 'ECS', 'ECA', 'TEC', 'EUU', 'FRO', 'PYF', 'GIB', ...
                  'GRL', 'GUM', 'HKG', 'IBD', 'ISL', 'IMN', 'IDA', 'IDX', 'IDA', ...
                  'IBT', 'IDB', 'LTE', 'LAC', 'TLA', 'LIC', 'LMY', 'LDC', 'MEA', ...
                  'MNA', 'TMN', 'MIC', 'OED', 'OSS', 'PSS', 'NAC', 'SXM', 'VIR', ...
                  'PSE', 'WLD', 'INX', 'PRE'};


Regions = {Canada, China, India, SouthKorea, Russia, Japan, US, EU, Africa, Brazil, ...
           OtherAsiaOceania, OtherEurope, LatinAmerica, Oceania, MiddleEast,OtherCountries};

pop_proj_region = cell(0, size(pop_proj, 2)); % Same number of columns as original data

for i = 1:length(Regions)
    region_countries = Regions{i};
    % Loop through the country codes of the current region
    for j = 1:length(region_countries)
        % Find the row in pop_proj corresponding to the current country code
        country_idx = find(strcmp(pop_proj(:,1), region_countries(j)));
        % Append the row to the new pop_proj_region array
        if ~isempty(country_idx)
            pop_proj_region = [pop_proj_region; pop_proj(country_idx, :)];
        end
    end
end

% Now, pop_proj_region will contain the population data reordered by the regions you defined.

pop_proj_region(:, [1]) = []

pop_proj=cell2mat(pop_proj_region)

% Define the region indices
region_counts = [1, 1, 2, 1, 1, 1, 1, 39, 51, 1, 47, 7, 32, 40, 49];

% Sum the population projections across each region
pop_proj_region_sum = cellfun(@(x) sum(x, 2), pop_proj_region(1:length(region_counts), :), 'UniformOutput', false);
pop_proj_region_sum=cell2mat(pop_proj_region_sum)
pop_proj_region_sum([15], :) = []

%%
var_y=CO2_Wfootprint_regionsTimeseries

Region14 = {"Canada", "China", "India", "SouthKorea", "Russia", "Japan", "US", "EU", "Africa", "Brazil", ...
           "OtherAsiaOceania", "OtherEurope", "LatinAmerica", "MiddleEast"}';

% Assuming CO2_Wfootprint_regionsTimeseries is a matrix of size 41x28
var_y1 = mat2cell(var_y, ones(49,1), ones(28,1));


Canada = [32];                % Canada
China = [31,41];                 % China
India = [35];                 % India
SouthKorea = [33];            % South Korea
Russia = [37];                % Russia
Japan = [30];                 % Japan
US = [29];                    % United States
EU = [1:28];                  % Austria, Belgium, Bulgaria, Cyprus, Czech Republic, Germany, Denmark, Estonia, Spain, Finland, France, Greece, Croatia, Hungary, Ireland, Italy, Lithuania, Luxembourg, Latvia, Malta, Netherlands, Poland, Portugal, Romania, Sweden, Slovenia, Slovakia, United Kingdom
Africa = [44, 48];            % South Africa, RoW Africa
Brazil = [34];                % Brazil
OtherAsiaOceania = [38,43,45];         % Australia, Indonesia
OtherEurope = [39, 40, 42, 47];  % Switzerland, Turkey, Norway, RoW Europe
LatinAmerica = [36, 46];      % Mexico, RoW America
MiddleEast = [49];            % RoW Middle East

regionOrder=[Canada, China, India, SouthKorea, Russia, Japan, US, EU, Africa, Brazil, ...
           OtherAsiaOceania, OtherEurope, LatinAmerica, MiddleEast]';
aggregationRegion=[1,2,1,1,1,1,1,28,2,1,3,4,2,1];
NNNR = 14; %New number of regions


m=1;
var_y2 = cell(size(var_y1));
for k=1:49
    for i=1:49
        if regionOrder(k,1)==i
            var_y2(m,:)=var_y1(i,:);
            m=m+1;
        end
    end
end

var_y3=mat2cell(var_y2,aggregationRegion,ones(28,1))
var_y4 = cellfun(@(x) sum(cell2mat(cellfun(@(y) sum(y(:), 'omitnan') * isnumeric(y), x, 'UniformOutput', false))), var_y3, 'UniformOutput', false);

var_y5=cell2mat(var_y4)
var_y6=[var_y5,zeros(14,28)]
var_y7=reshape(var_y6,[],1)


%% Yearid and Regionid
year=reshape(repmat(1995:1:2050,14,1),[],1)
yearid=categorical(reshape(repmat(1:1:56,14,1),[],1));

regionid_1=repmat(Region14,1,56)
regionid_2=reshape(regionid_1,[],1)

regionid_3=1:1:14
regionid_4=repmat(regionid_3,1,56)
regionid_5=reshape(regionid_4,[],1)
% establish table
regionid=categorical(regionid_5)
regionname=regionid_2

id=categorical([1:784]');

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
Region14 = {"Canada", "China", "India", "SouthKorea", "Russia", "Japan", "US", "EU", "Africa", "Brazil", ...
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
    title(Region14{i}, 'FontSize', 8);  % Use the corresponding region name from Region14
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

% Get the list of unique regions
unique_regions = unique(RegTable.regionid);

% Define the region names in the correct order
Region14_title = {"Canada", "China", "India", "SouthKorea", "Russia", "Japan", "US", "EU", "Africa", "Brazil", ...
           "Other Asia&Oceania", "Other Europe", "Latin America", "Middle East"}';

% Convert yearid to actual years (1995-2050)
% Assume that yearid starts at 1 for 1995 and ends at 56 for 2050
years = linspace(1995, 2050, 56);

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
    fill([years(29:end)'; flipud(years(29:end)')], ...
         [future_data.lower_bound; flipud(future_data.upper_bound)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
     
    % Plot the predicted carbon footprint (historical + future)
    plot(years, region_data.predicted_carbon, 'r-', 'LineWidth', 2);
    
    % Plot the actual historical carbon footprint
    plot(years, region_data.y, 'k--', 'LineWidth', 2);
    
    % Set title and labels
    title(Region14_title{i}, 'FontSize', 14);  % Use the corresponding region name from Region14
    
    % Standardize y-tick format to 10^11 for all subplots
    
    % Set custom ylim for Africa and US
    if strcmp(Region14_title{i}, 'Africa')
        ylim([0 7e11]);  
    elseif strcmp(Region14_title{i}, 'Japan')
        ylim([0 5e11]);  
    elseif strcmp(Region14_title{i}, 'SouthKorea')
        ylim([0 3e11]);
    elseif strcmp(Region14_title{i}, 'EU')
        ylim([0 2e12]); 
    elseif strcmp(Region14_title{i}, 'Brazil')
        ylim([0 4e11]);  
    elseif strcmp(Region14_title{i}, 'US')
        ylim([3e11 12e11]);  % Custom ylim for US
    end

    xlim([1995 2050])

    hold off;
end

set(gcf, 'Position',  [225,99,1127,902])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/construction/figures/fullregion_projection.svg');

