Code package for paper:
Chaohui Li, Prajal Pradhan∗, Guoqian Chen∗, Juergen Kropp, and Hans Joachim Schellnhuber. Carbon Cost of Constructing the Global Built Environment Can Preclude the Climate Change Goals. 

All data used in this paper are publicly available. Input-Output Tables used in this study are publicly available at https://www.exiobase.eu and https://zenodo.org/records/5589597, socio-economic data of historical and future projections are available at https://iiasa.ac.at/models-tools-data/ssp,  https://data.worldbank.org, latest carbon budget data are available at https://essd.copernicus.org/articles/15/2295/2023/essd-15-2295-2023.html 

construction_ini.m is a file that contains the initial settings. This code should be run prior to all other codes. This code defines the sector and region aggregation in the input-output analysis. It also defines the developed, developing, and emerging economies.

construction_footprints.m conducts input-output analysis and calculates the historical trend of the construction sector from 1995 to 2022. The output of this script is the global evolution trend of the construction industry, with detailed sectoral information, supply-chain information, region-specific information, and other details.

construction_regression_SI.m should be run prior to series of scripts for future projections. This script conducts a series of statistical tests in regression to avoid multicollinearity, autocorrelation, heteroscedasticity, etc. This code file includes tests such as KPSS (Kwiatkowski-Phillips-Schmidt-Shin), PP tests, Augmented Dickey-Fuller (ADF), Pearson Correlation Test, Ljung-Box Test, Durbin-Watson Test, ARCH and GARCH test.

construction_regression_country_carbon.m calculates the country-specific regressions for all future SSPs scenarios.

construction_projection_fullcountry9.m divides the world economy into 9 regions and projects future trajectories of 9 regions for the business-as-usual scenario.

construction_prediction_fullcountry_fixedeffect.R and construction_prediction_fullcountry.R uses OLS, fixed effect and ARIMA to project future trajectories.

construction_carbon_budget_2degrees.m and construction_carbon_budget15degrees.m uses the exponential decay model to calculate the carbon budget trajectories and juxtapose with the projected future pathways of the construction industry to map intersection points.
