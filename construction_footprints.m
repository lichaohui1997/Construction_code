%% load data
for yearnum=1995:2022;
    year=num2str(yearnum);
% data repository
file=strcat('IOT_',year,'_ixi');
addpath(fullfile('exiodata','ixi',file));
addpath(fullfile('exiodata','ixi',file,'satellite'));
% load data
A_=readmatrix('A.txt');%A=coefficient matrix of intermediate matrix(EXIOBASE)
A=A_(2:7988,3:7989);
F_=readmatrix('Y.txt');%Y=final demand matrix of raw data(EXIOBASE), in this script denoted as "F"
F=F_(2:7988,3:345);
S_=readmatrix('F.txt');%F is the satellite account matrix of raw data(EXIOBASE), in this script denoted as "S"
S=S_(:,2:7988);
P=S(1:9,:);%P=primary input matrix
% Use coeficient A to calculate intermediate matrix T (for EXIOBASE dataset only)
Fsum2=sum(F,2);
I=eye(NI*NR);
B=(I-A);
X1=(B'*B)\(B'*Fsum2);%equivalent to pinv(B) 
diagX1=diag(X1);
T=A*diagX1;
P(find(isnan(P)==1)) = 0
% load data
DATA=[24,93,94,428,438,439];%Line numbers of carbon emissions in satelite account of IO table
D_=S(DATA,:);
D=sum(D_,1)';
%% Solve for intensity
X2=sum(T,2)+sum(F,2);
diagX2=diag(X2);
Prow=sum(P);
Psum=sum(sum(P),2);
Household=zeros(NI*NR,NR);
Government=zeros(NI*NR,NR);
Non_profit=zeros(NI*NR,NR);
for i=1:NR
    Household(:,i)=F(:,(i-1)*NF+1);
    Non_profit(1:NI*NR,i)=F(1:NI*NR,(i-1)*NF+2);
    Government(1:NI*NR,i)=F(1:NI*NR,(i-1)*NF+3); 
end
Consumption=Household+Non_profit+Government;
Frest=sum(F,2)-sum(Consumption,2);
FP=Frest*Prow;
FP1=FP/Psum;
C=diagX2-T-FP1;
intensity=D'*pinv(C);
EECsum=sum(sum(Consumption,2).*intensity');
DEEsum=sum(D',2);%IF EECsum=DEEsum, then the calculation up to here is correct.
%% Resources embodied in final demand
Fembodiedenergy1=F.*repmat(intensity',1,NR*NF);
Fembodiedenergy2=mat2cell(Fembodiedenergy1,ones(NR*NI,1)*1,ones(NR,1)*NF)
Fembodiedenergy=cell2mat(cellfun(@sum,Fembodiedenergy2,'UniformOutput',false));
%% Resources embodied in consumption
Consumptionembodiedenergy=Consumption.*repmat(intensity',1,NR);
Consumptionembodiedenergysum=sum(Consumptionembodiedenergy,2);
Consumptionembodiedenergysector=sum(reshape(Consumptionembodiedenergysum,NI,NR),2);
%% Resources embodied in intermediate trade
Tembodiedenergy=T.*repmat(intensity',1,NI*NR);
%% Aggregation of sectors in input-output tables, intermediate matrix
hh1=mat2cell(Tembodiedenergy,ones(NR,1)*NI,ones(NR,1)*NI);
hh2=cat(3,hh1{:});
% sort out fist dimension
m=1;
hh3=zeros(NI,NI,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            hh3(m,:,:)=hh2(i,:,:);
            m=m+1;
        end
    end
end
% sort out second dimension
m=1;
hh4=zeros(NI,NI,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            hh4(:,m,:)=hh3(:,i,:);
            m=m+1;
        end
    end
end
% aggregate sectors
hh5=mat2cell(hh4,aggregation,aggregation,ones(NR*NR/1,1)*1);
hh6=cellfun(@sum,cellfun(@sum,hh5,'UniformOutput',false));
hh7=inverse_cat(1,hh6);
hh8=mat2cell(hh7,ones(1/1,1)*NNI,ones(NR*NR,1)*NNI);
hh9=reshape(hh8,NR,NR);
% sort out first dimension 
m=1;
hh10= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            hh10(m,:)=hh9(i,:);
            m=m+1;
        end
    end
end
% sort out second dimension
m=1;
hh11= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            hh11(:,m)=hh10(:,i);
            m=m+1;
        end
    end
end
% aggregate sectors
hh12=mat2cell(hh11,aggregationCountry,aggregationCountry)
hh13=cell(NNR,NNR)

for i=1:NNR
    for k=1:NNR
        hh13{i,k}=sum(cat(3,hh12{i,k}{:}),3);
        if i==k
            hh13{i,k}=zeros(NNI,NNI)
        end
    end
end
for i=1:NNR
    for k=1:NNR
        hh14{i,k}=sum(cat(3,hh12{i,k}{:}),3);
       
    end
end
hh15=cell2mat(hh14)
%% Aggregation of sectors in input-output tables in final demand matrix
% final demand
jj1=mat2cell(Fembodiedenergy,ones(NR,1)*NI,ones(NR,1)*1);
jj2=cat(3,jj1{:});
%
m=1;
jj3=zeros(NI,1,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            jj3(m,:,:)=jj2(i,:,:);
            m=m+1;
        end
    end
end

jj4=mat2cell(jj3,aggregation,1,ones(NR*NR/1,1)*1);
jj5=cellfun(@sum,cellfun(@sum,jj4,'UniformOutput',false));
jj6=inverse_cat(1,jj5);
jj7=reshape(jj6,NR*NNI,NR)
jj8=mat2cell(jj7,ones(NR*NNI/NNI,1)*NNI,ones(NR/1,1)*1);
%
m=1;
j9= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            j9(m,:)=jj8(i,:);
            m=m+1;
        end
    end
end
%
m=1;
jj10= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            jj10(:,m)=j9(:,i);
            m=m+1;
        end
    end
end
%
jj11=mat2cell(jj10,aggregationCountry,aggregationCountry)
jj12=cell(NNR,NNR)
for i=1:NNR
    for k=1:NNR
        jj12{i,k}=sum(cat(3,jj11{i,k}{:}),3);
        if i==k
            jj12{i,k}=zeros(NNI,1)
        end
    end
end

for i=1:NNR
    for k=1:NNR
        jj13{i,k}=sum(cat(3,jj11{i,k}{:}),3);
     
    end
end
jj14=cell2mat(jj13)
miniZF=[hh15,jj14];
WminiZF=miniZF;
fig1_1=sum(miniZF,2);
fig1=reshape(fig1_1,[NNI,NNR]);

%% This section of codes outputs the circos diagram
% Embodied emissions/energyuse in intermediate matrix
h1=mat2cell(Tembodiedenergy,ones(NR,1)*NI,ones(NR,1)*NI);
h2=cat(3,h1{:});
m=1;
h3=zeros(NI,NI,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            h3(m,:,:)=h2(i,:,:);
            m=m+1;
        end
    end
end
m=1;
h4=zeros(NI,NI,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            h4(:,m,:)=h3(:,i,:);
            m=m+1;
        end
    end
end
h5=mat2cell(h4,aggregation,aggregation,ones(NR*NR/1,1)*1);
h6=cellfun(@sum,cellfun(@sum,h5,'UniformOutput',false));
h6(NNI,:,:)=[];
h6(:,NNI,:)=[];
h7=inverse_cat(1,h6);
h8=mat2cell(h7,ones(1/1,1)*NNI-1,ones(NR*NR*(NNI-1)/(NNI-1),1)*(NNI-1));
h9=reshape(h8,NR,NR);

m=1;
h10= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            h10(m,:)=h9(i,:);
            m=m+1;
        end
    end
end
%
m=1;
h11= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            h11(:,m)=h10(:,i);
            m=m+1;
        end
    end
end

h12=mat2cell(h11,aggregationCountry,aggregationCountry)
h13=cell(NNR,NNR)
for i=1:NNR
    for k=1:NNR
        h13{i,k}=sum(cat(3,h12{i,k}{:}),3);
        if i==k
            h13{i,k}=zeros(NNI-1,NNI-1)
        end
    end
end
for i=1:NNR
    for k=1:NNR
        h14{i,k}=sum(cat(3,h12{i,k}{:}),3);
       
    end
end
h_firstmatrix=cellfun(@sum,cellfun(@sum,h13,'UniformOutput',false));
h_firstmatrix_withselftrade=cellfun(@sum,cellfun(@sum,h14,'UniformOutput',false));
%% % Embodied emissions/energyuse in Final demand matrix
j1=mat2cell(Fembodiedenergy,ones(NR,1)*NI,ones(NR,1)*1);
j2=cat(3,j1{:});
m=1;
j3=zeros(NI,1,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            j3(m,:,:)=j2(i,:,:);
            m=m+1;
        end
    end
end

j4=mat2cell(j3,aggregation,1,ones(NR*NR/1,1)*1);
j5=cellfun(@sum,cellfun(@sum,j4,'UniformOutput',false));
j5(NNI,:,:)=[];
j6=inverse_cat(1,j5);
j7=reshape(j6,NR*(NNI-1),NR)
j8=mat2cell(j7,ones(NR*(NNI-1)/(NNI-1),1)*(NNI-1),ones(NR/1,1)*1);

m=1;
j9= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            j9(m,:)=j8(i,:);
            m=m+1;
        end
    end
end
%
m=1;
j10= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if countryOrder(k,1)==i
            j10(:,m)=j9(:,i);
            m=m+1;
        end
    end
end
%
j11=mat2cell(j10,aggregationCountry,aggregationCountry)
j12=cell(NNR,NNR)
for i=1:NNR
    for k=1:NNR
        j12{i,k}=sum(cat(3,j11{i,k}{:}),3);
        if i==k
            j12{i,k}=zeros(NNI-1,1)
        end
    end
end

for i=1:NNR
    for k=1:NNR
        j13{i,k}=sum(cat(3,j11{i,k}{:}),3);
     
    end
end
j_secondmatrix=cellfun(@sum,cellfun(@sum,j12,'UniformOutput',false));
j_secondmatrix_withselftrade=cellfun(@sum,cellfun(@sum,j13,'UniformOutput',false));
regionIO=h_firstmatrix+j_secondmatrix;
regionIO_withselftrade=h_firstmatrix_withselftrade+j_secondmatrix_withselftrade;
%% This section calculates the embodied emissions in intermediate matrix
h1=mat2cell(Tembodiedenergy,ones(NR,1)*NI,ones(NR,1)*NI);
h2=cat(3,h1{:});
% sectoral aggregation 
m=1;
h3=zeros(NI,NI,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            h3(m,:,:)=h2(i,:,:);
            m=m+1;
        end
    end
end

m=1;
h4=zeros(NI,NI,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            h4(:,m,:)=h3(:,i,:);
            m=m+1;
        end
    end
end

h5=mat2cell(h4,aggregation,aggregation,ones(NR*NR/1,1)*1);
h6=cellfun(@sum,cellfun(@sum,h5,'UniformOutput',false));
h6(NNI,:,:)=[];
h6(:,NNI,:)=[];
h7=inverse_cat(1,h6);
h8=mat2cell(h7,ones(1/1,1)*NNI-1,ones(NR*NR*(NNI-1)/(NNI-1),1)*(NNI-1));
h9=reshape(h8,NR,NR);

m=1;
h10= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if regionOrder(k,1)==i
            h10(m,:)=h9(i,:);
            m=m+1;
        end
    end
end
%
m=1;
h11= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if regionOrder(k,1)==i
            h11(:,m)=h10(:,i);
            m=m+1;
        end
    end
end
%
h12=mat2cell(h11,aggregationRegion,aggregationRegion)
%
h13=cell(NNNR,NNNR)

for i=1:NNNR
    for k=1:NNNR
        h13{i,k}=sum(cat(3,h12{i,k}{:}),3);
        if i==k
            h13{i,k}=zeros(NNI-1,NNI-1)
        end
    end
end

for i=1:NNNR
    for k=1:NNNR
        h14{i,k}=sum(cat(3,h12{i,k}{:}),3);
       
    end
end
%
h_firstmatrix=cellfun(@sum,cellfun(@sum,h13,'UniformOutput',false));
h_firstmatrix_withselftrade=cellfun(@sum,cellfun(@sum,h14,'UniformOutput',false));
%% This section of codes calculates the embodied emissions in the final demand matrix
j1=mat2cell(Fembodiedenergy,ones(NR,1)*NI,ones(NR,1)*1);
j2=cat(3,j1{:});

%
m=1;
j3=zeros(NI,1,NR*NR);
for k=1:NI
    for i=1:NI
        if Order(k,1)==i
            j3(m,:,:)=j2(i,:,:);
            m=m+1;
        end
    end
end

j4=mat2cell(j3,aggregation,1,ones(NR*NR/1,1)*1);
j5=cellfun(@sum,cellfun(@sum,j4,'UniformOutput',false));
j5(NNI,:,:)=[];
j6=inverse_cat(1,j5);
j7=reshape(j6,NR*(NNI-1),NR)
j8=mat2cell(j7,ones(NR*(NNI-1)/(NNI-1),1)*(NNI-1),ones(NR/1,1)*1);

m=1;
j9= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if regionOrder(k,1)==i
            j9(m,:)=j8(i,:);
            m=m+1;
        end
    end
end

m=1;
j10= cell(NR,NR)
for k=1:NR
    for i=1:NR
        if regionOrder(k,1)==i
            j10(:,m)=j9(:,i);
            m=m+1;
        end
    end
end
%
j11=mat2cell(j10,aggregationRegion,aggregationRegion)
j12=cell(NNNR,NNNR)

for i=1:NNNR
    for k=1:NNNR
        j12{i,k}=sum(cat(3,j11{i,k}{:}),3);
        if i==k
            j12{i,k}=zeros(NNI-1,1)
        end
    end
end

for i=1:NNNR
    for k=1:NNNR
        j13{i,k}=sum(cat(3,j11{i,k}{:}),3);
     
    end
end
%
j_secondmatrix=cellfun(@sum,cellfun(@sum,j12,'UniformOutput',false));
j_secondmatrix_withselftrade=cellfun(@sum,cellfun(@sum,j13,'UniformOutput',false));
regionIO_8=h_firstmatrix+j_secondmatrix;
regionIO_withselftrade_8=h_firstmatrix_withselftrade+j_secondmatrix_withselftrade;

%%
p1=mat2cell(hh15,ones(NNR,1)*NNI,ones(NNR,1)*NNI)
p2=cat(3,p1{:});
p3=p2(1:NNI,N_Construction,:)
Sector_all=sum(p3(:,:),3);
Sector_1=mat2cell(Sector_all,ones(1,1)*13,ones(NR,1)*NR)
Sector_2=cell2mat(cellfun(@(x) sum(x, 2),Sector_1,'UniformOutput',false));
Sector_sum=sum(Sector_2)
Sector_primary=fig1(12,:)-Sector_sum;
%% This calculates the capital investments
Sector_3=[Sector_2(1:12,:);zeros(1,49);Sector_2(13,:)] 
Sector_3(13,:)=Sector_primary
Sector_4=Sector_3./sum(Sector_3) 
%% These are the output information 
Wfig1=fig1;
Wfootprint_construction=sum(Wfig1(12,:),2);
Wprct_construction=Wfootprint_construction/sum(D);
Wprct_construction=Wfootprint_construction/sum(D);
Wfootprint_regions=Wfig1(12,:);
Wp1=p1;
WsumD=sum(D);
Wsector=Sector_3';
Wsector_prct=Sector_4';
Wfootprint_types=sum(Wsector)';
Wfootprint_types_prct=Wfootprint_types./sum(Wfootprint_types)
WminiIO_Newregion=regionIO;
WregionImport=sum(WminiIO_Newregion);
WregionExport=-1*sum(WminiIO_Newregion,2)';
WImEx=[WregionExport;WregionImport]';
WImEx22=[sum(WImEx(1:28,:));WImEx(29:49,:)];
WBalance=[WregionImport+WregionExport]';
WBalance22=[sum(WBalance(1:28,:));WBalance(29:49,:)];
Wcircos=regionIO_8;
WregionImport=sum(WminiIO_Newregion);
WregionExport=-1*sum(WminiIO_Newregion,2)';
WImEx=[WregionExport;WregionImport]';
WImEx22=[sum(WImEx(1:28,:));WImEx(29:49,:)];
Wintensity=intensity_7';
Wintensity_all=intensity_8';
end
 