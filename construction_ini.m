tic
clc;clear;close all
% by Chaohui Li
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/construction'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/other'));

%% Initializations

NI = 163;%Number of industries
NR = 49;%Number of regions
NF = 7;%Number of final demand sectors

[~,regioncodes,~] = xlsread('EXIOBASE_metadata.xlsx','Countries','B1:B49');
[~,regionnames,~] = xlsread('EXIOBASE_metadata.xlsx','Countries','C1:C49');
[ProdIndConcordance,~,~] = xlsread('EXIOBASE_metadata.xlsx','ProdIndConcordance','F7:FL206');
[population,~,~] = xlsread('population.xlsx','data','A1:AJ8');
ssp_prct=xlsread('ssp_emissions.xlsx','data','G14:AQ18');

%% Industry aggregation
% Define your specific categories
Steel = [25, 72, 73]; %3
Cement = [69];%1
Glass = [65, 66];%2
OtherMetal = [24, 26:31, 74:85];%1+6+12=19
Transport = [120:126];%7
CeramicClay = [67, 68];%2
Service = [115:119, 127:138, 159:163];%5+12+5=22
CPR = [63:64];% ChemicalPlasticRubber %2
ClinkerAsh = [70,71];%2
Machinery = [86:93];%8
Biobased = [7, 8, 15, 18, 50, 54];%6
Construction = [113, 114];%2

combinedCategories = [Steel, Cement, Glass, OtherMetal, Transport, CeramicClay, Service, CPR, ClinkerAsh, Machinery, Biobased, Construction];
fullRange = 1:163;
Industry = setdiff(fullRange, combinedCategories);

Order=[Cement, ClinkerAsh,CeramicClay,Steel,OtherMetal, Glass, CPR,  Biobased, Transport, Service, Machinery, Construction, Industry]';
aggregation=[1,2,2,3,19,2,2,6,7,22,8,2,87];

NNI = 13; %New number of industries

N_Cement=1;
N_ClinkerAsh=2;
N_CeramicClay=3;
N_Steel=4;
N_OtherMetal=5;
N_Glass=6;
N_CPR=7;
N_Biobased=8;
N_Transport=9;
N_Service=10;
N_Machinery=11;
N_Construction=12;
N_Industry=13;
%% Region aggregation 49
%%
countryOrder= readmatrix('countryorder.xlsx','Sheet',1,'Range','A:A');
aggregationCountry= readmatrix('countryorder.xlsx','Sheet',2,'Range','A:A')';
countryNameFull=readvars('countryorder.xlsx','Sheet',1,'Range','C1:C49');
countryName= readvars('countryorder.xlsx','Sheet',1,'Range','E1:E49');
countryNameMain= readvars('countryorder.xlsx','Sheet',1,'Range','F1:F49');
countryName22=readvars('countryorder.xlsx','Sheet',3,'Range','A:A');
GDP= readmatrix('population.xlsx','Sheet','GDP','Range','B2:B50');
NNR = 49; %New number of regions
%%
%% Region aggregation 8
EU=[1:28];
US=[29];
China=[31,41];
OtherAsiaPacific=[30,33,35,38,43,45];
MiddleEast=[49];
Africa=[44,48];
OtherEurope=[37,39,40,42,47];
OtherAmerica=[32,34,36,46];
regionOrder=[EU,US,China,OtherAsiaPacific,MiddleEast,Africa,OtherEurope,OtherAmerica]';
aggregationRegion=[28,1,2,6,1,2,5,4];
NNNR = 8; %New number of regions

N_EU=1;
N_US=2;
N_China=3;
N_OtherAsiaPacific=4;
N_MiddleEast=5;
N_Africa=6;
N_OtherEurope=7;
N_OtherAmerica=8;
RegionNames={'EU','US','China','Other Asia&Pacific','Middle East','Africa','Other Europe&Eurasia','Other America'}

developed=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,33,38,39,41,42];
emerging=[31,34,35,36,37,40,43,44];
lowincome=[45,46,47,48,49];
developOrder=[developed,emerging,lowincome];
aggregationDevelop=[36,8,5];

developOrder2=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,2,2,2,1,1,2,2,1,2,2,3,3,3,3,3]';

%%
RegionNames={'EU','US','China','Other Asia&Pacific','Middle East','Africa','Other Europe&Eurasia','Other America'}
constructionNames={'Construction'}
