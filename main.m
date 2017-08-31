%30 June 2017
%Taisa Kushner
%fits for data

%turn off the beeping
beep off 
clear all

%pull all files for given patient
pID='PSO3-001-0001';

%% import food data to clean up files
[~, ~, raw] = xlsread('/Users/tkushner/Documents/AP_Controller/PSO3Data/PSO3Food.xlsx','Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[2,5,6]);
raw = raw(:,[1,3,4]);

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
Food.SessionID = data(:,1);
Food.PtID = cellVectors(:,1);
Food.MealDtTm1 = data(:,2);
Food.GramsCHO1 = data(:,3);
Food.ForLowTrt1 = cellVectors(:,2);
Food.ProteinAndFat1 = cellVectors(:,3);

%% pull session numbers that have bolus for given patient
Pull=strncmp(Food.PtID,pID,13);
sesToPull=unique(Food.SessionID(Pull));

%%
allfiles=dir(strcat('../outputs/byPatient/session-',pID,'-*'));
%how many are there
numfiles=size(allfiles);
MAX=numfiles(1);

%Initialize patient struct for speed
patient(MAX).SessionID=[];
patient(MAX).Datetime=[];
patient(MAX).CGM=[];
patient(MAX).IOB=[];
patient(MAX).Bolus=[];
patient(MAX).BkgInsulin=[];
patient(MAX).times=[];
patient(MAX).gtimes=[];
patient(MAX).gCGM=[];
patient(MAX).gIOB=[];
patient(MAX).gDatetime=cell(1,1);

tic
parfor i=1:MAX
    [patient(i).SessionID,patient(i).Datetime,patient(i).CGM,patient(i).IOB,patient(i).Bolus,patient(i).BkgInsulin] = ...
        importCGMDATA(strcat('../outputs/byPatient/',allfiles(i).name));
    
    formatin='HH:MM:SS mm/dd/yyyy';
    patient(i).times=datenum(patient(i).Datetime,formatin);
    %rescale times based on start time
    patient(i).times=patient(i).times-patient(i).times(1);
    %get timestep - keep 5min intervals (5min~.0035)
    step=diff(patient(i).times);
    patient(i).times=patient(i).times./step(2)*5;
    
    %pull only times after first 180min burn-in
    gindex=find((patient(i).times>=180)-ismissing(patient(i).IOB));
    patient(i).gtimes=round(patient(i).times(gindex));
    patient(i).gCGM=patient(i).CGM(gindex);
    patient(i).gIOB=patient(i).IOB(gindex);
    patient(i).gDatetime=patient(i).Datetime(gindex);
    
    if isempty(patient(i).gtimes)
        patient(i).gtimes=[0,0];
        patient(i).gCGM=0;
        patient(i).gIOB=0;
    end
end
toc

%% clean data
n=1;
for i=1:MAX
    if ~ismember(patient(i).SessionID(1),sesToPull)
        cleanpatient(n)=patient(i);
        n=n+1;
    end
end

%% Optimization by window
tic 
MAX=numel(cleanpatient);
a0=[.6 .4 -2.5];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=6;
Delta1=6;
stepsz=24;
ovlp=5;

[modelFits30Win24, stats30Win24] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);
%[modelFits30Win24_dirty, stats30Win24_dirty]=WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

stepsz=36;
ovlp=5;
[modelFits30Win36, stats30Win36] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);

stepsz=60;
ovlp=5;
[modelFits30Win60, stats30Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);

% 45min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=9;
Delta1=9;

stepsz=24;
ovlp=5;
[modelFits45Win24, stats45Win24] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);


stepsz=36;
ovlp=5;
[modelFits45Win36, stats45Win36] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits45Win60, stats45Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);

% 60min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=12;
Delta1=12;

stepsz=36;
ovlp=5;
[modelFits60Win36, stats60Win36] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits60Win60, stats60Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);

%120min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=24;
Delta1=24;

stepsz=60;
ovlp=5;

[modelFits120Win60, stats120Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX,stepsz,ovlp);
toc
%% plots for presentation

Delta1=6;
Delta2=9;
Delta3=12;
Delta4=24;

PlotsRange=23:27;
%%
for i=PlotsRange
    A1=[patient(i).gCGM(1:end-2*Delta1), patient(i).gCGM(Delta1+1:end-Delta1), patient(i).gIOB(1:end-2*Delta1)];
    A2=[patient(i).gCGM(1:end-2*Delta2), patient(i).gCGM(Delta2+1:end-Delta2), patient(i).gIOB(1:end-2*Delta2)];
    
    %x=[modelFits30Win24(i).MEAN(1); modelFits30Win24(i).MEAN(2); modelFits30Win24(i).MEAN(3)];
    x1=[stats30Win24.mean(1); stats30Win24.mean(2); stats30Win24.mean(3)];
    x2=[stats45Win24.mean(1); stats45Win24.mean(2); stats45Win24.mean(3)];
    predict1=A1*x1;
    predict2=A2*x2;
    predictmin1=predict1-stats30Win24.RES95;
    predictmax1=predict1+stats30Win24.RES95;
    
    predictmin2=predict2-stats45Win24.RES95;
    predictmax2=predict2+stats45Win24.RES95;
    
    predictMAX=zeros(size(predict1));
    predictMIN=zeros(size(predict1));
    diffdelta=2*Delta2-2*Delta1;
    for n=1:diffdelta
        predictMAX(n)=predictmax1(n);
        predictMIN(n)=predictmin1(n);
    end
    
    for n=diffdelta+1:length(predictMAX)
        predictMAX(n)=min(predictmax1(n),predictmax2(n-diffdelta));
        predictMIN(n)=max(predictmin1(n),predictmin2(n-diffdelta));
    end
    
    
    figure(i)
    subplot(2,1,1)
    %plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gCGM(2*Delta1+1:end),...
    %    patient(i).gtimes(2*Delta1+1:end),predictmin1,'r--',patient(i).gtimes(2*Delta1+1:end),predictmax1,'r--')
    hold on
    %plot(patient(i).gtimes(2*Delta2+1:end),predictmin2,'r--',patient(i).gtimes(2*Delta2+1:end),predictmax2,'r--')
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gCGM(2*Delta1+1:end),...
        patient(i).gtimes(2*Delta1+1:end),predictMIN,'r--',patient(i).gtimes(2*Delta1+1:end),predictMAX,'r--')
    xlabel('time since start trial','FontSize',14)
    ylabel('Glucose (mg/dL)','FontSize',14)
    legend({'CGM reading','Predicted range, 99% confidence'},'FontSize',14)
    title(strcat('120 min Window PSO3-001-0001 trial',num2str(i)),'FontSize',16)
    subplot(2,1,2)
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gIOB(2*Delta1+1:end))
    xlabel('time since start trial','FontSize',14)
    ylabel('IOB (units)','FontSize',14)
end
%%

for i=PlotsRange
    A1=[patient(i).gCGM(1:end-2*Delta1), patient(i).gCGM(Delta1+1:end-Delta1), patient(i).gIOB(1:end-2*Delta1)];
    A2=[patient(i).gCGM(1:end-2*Delta2), patient(i).gCGM(Delta2+1:end-Delta2), patient(i).gIOB(1:end-2*Delta2)];
    A3=[patient(i).gCGM(1:end-2*Delta3), patient(i).gCGM(Delta3+1:end-Delta3), patient(i).gIOB(1:end-2*Delta3)];

    
    x1=[stats30Win36.mean(1); stats30Win36.mean(2); stats30Win36.mean(3)];
    x2=[stats45Win36.mean(1); stats45Win36.mean(2); stats45Win36.mean(3)];
    x3=[stats60Win36.mean(1); stats60Win36.mean(2); stats60Win36.mean(3)];
    predict1=A1*x1;
    predict2=A2*x2;
    predict3=A3*x3;
    predictmin1=predict1-stats30Win36.RES95;
    predictmax1=predict1+stats30Win36.RES95;
    
    predictmin2=predict2-stats45Win36.RES95;
    predictmax2=predict2+stats45Win36.RES95;
    
    predictmin3=predict3-stats60Win36.RES95;
    predictmax3=predict3+stats60Win36.RES95;
    
    predictMAX=zeros(size(predict1));
    predictMIN=zeros(size(predict1));
    diffdelta=2*Delta2-2*Delta1;
    diffdelta2=2*Delta3-2*Delta1;
    
    for n=1:diffdelta
        predictMAX(n)=predictmax1(n);
        predictMIN(n)=predictmin1(n);
    end
        
    for n=diffdelta+1:diffdelta2
        predictMAX(n)=min(predictmax1(n),predictmax2(n-diffdelta));
        predictMIN(n)=max(predictmin1(n),predictmin2(n-diffdelta));
    end
    
    for n=diffdelta2+1:length(predictMAX)
        predictMAX(n)=min([predictmax1(n),predictmax2(n-diffdelta),predictmax3(n-diffdelta2)]);
        predictMIN(n)=max([predictmin1(n),predictmin2(n-diffdelta),predictmin3(n-diffdelta2)]);
    end
    
    figure(100+i)
    subplot(2,1,1)
    %plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gCGM(2*Delta1+1:end),...
    %    patient(i).gtimes(2*Delta1+1:end),predictmin1,'r--',patient(i).gtimes(2*Delta1+1:end),predictmax1,'r--')
    hold on
    %plot(patient(i).gtimes(2*Delta2+1:end),predictmin2,'r--',patient(i).gtimes(2*Delta2+1:end),predictmax2,'r--')
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gCGM(2*Delta1+1:end),...
        patient(i).gtimes(2*Delta1+1:end),predictMIN,'r--',patient(i).gtimes(2*Delta1+1:end),predictMAX,'r--')
    xlabel('time since start trial','FontSize',14)
    ylabel('Glucose (mg/dL)','FontSize',14)
    legend({'CGM reading','Predicted range, 99% confidence'},'FontSize',14)
    title(strcat('180min Window PSO3-001-0001 trial',num2str(i)),'FontSize',16)
    subplot(2,1,2)
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gIOB(2*Delta1+1:end))
    xlabel('time since start trial','FontSize',14)
    ylabel('IOB (units)','FontSize',14)
end

%%
for i=PlotsRange
    A1=[patient(i).gCGM(1:end-2*Delta1), patient(i).gCGM(Delta1+1:end-Delta1), patient(i).gIOB(1:end-2*Delta1)];
    A2=[patient(i).gCGM(1:end-2*Delta2), patient(i).gCGM(Delta2+1:end-Delta2), patient(i).gIOB(1:end-2*Delta2)];
    A3=[patient(i).gCGM(1:end-2*Delta3), patient(i).gCGM(Delta3+1:end-Delta3), patient(i).gIOB(1:end-2*Delta3)];

    x1=[stats30Win60.mean(1); stats30Win60.mean(2); stats30Win60.mean(3)];
    x2=[stats45Win60.mean(1); stats45Win60.mean(2); stats45Win60.mean(3)];
    x3=[stats60Win36.mean(1); stats60Win36.mean(2); stats60Win36.mean(3)];

    predict1=A1*x1;
    predict2=A2*x2;
    predict3=A3*x3;

    predictmin1=predict1-stats30Win60.RES95;
    predictmax1=predict1+stats30Win60.RES95;
    
    predictmin2=predict2-stats45Win60.RES95;
    predictmax2=predict2+stats45Win60.RES95;
    
    predictmin3=predict3-stats60Win36.RES95;
    predictmax3=predict3+stats60Win36.RES95;
    
    predictMAX=zeros(size(predict1));
    predictMIN=zeros(size(predict1));
    diffdelta=2*Delta2-2*Delta1;
    for n=1:diffdelta
        predictMAX(n)=predictmax1(n);
        predictMIN(n)=predictmin1(n);
    end
    
    for n=diffdelta+1:diffdelta2
        predictMAX(n)=min(predictmax1(n),predictmax2(n-diffdelta));
        predictMIN(n)=max(predictmin1(n),predictmin2(n-diffdelta));
    end
    
    for n=diffdelta2+1:length(predictMAX)
        predictMAX(n)=min([predictmax1(n),predictmax2(n-diffdelta),predictmax3(n-diffdelta2)]);
        predictMIN(n)=max([predictmin1(n),predictmin2(n-diffdelta),predictmin3(n-diffdelta2)]);
    end
    
    
    figure(200+i)
    subplot(2,1,1)
    %plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gCGM(2*Delta1+1:end),...
    %    patient(i).gtimes(2*Delta1+1:end),predictmin1,'r--',patient(i).gtimes(2*Delta1+1:end),predictmax1,'r--')
    hold on
    %plot(patient(i).gtimes(2*Delta2+1:end),predictmin2,'r--',patient(i).gtimes(2*Delta2+1:end),predictmax2,'r--')
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gCGM(2*Delta1+1:end),...
        patient(i).gtimes(2*Delta1+1:end),predictMIN,'r--',patient(i).gtimes(2*Delta1+1:end),predictMAX,'r--')
    xlabel('time since start trial','FontSize',14)
    ylabel('Glucose (mg/dL)','FontSize',14)
    legend({'CGM reading','Predicted range, 99% confidence'},'FontSize',14)
    title(strcat('300min Window PSO3-001-0001 trial',num2str(i)),'FontSize',16)
    subplot(2,1,2)
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gIOB(2*Delta1+1:end))
    xlabel('time since start trial','FontSize',14)
    ylabel('IOB (units)','FontSize',14)
end

%% write to txt for tex
% writetoFileTex('30Win24',pID,stats30Win24);
% writetoFileTex('30Win36',pID,stats30Win36);
% writetoFileTex('30Win60',pID,stats30Win60);
% writetoFileTex('45Win24',pID,stats45Win24);
% writetoFileTex('45Win36',pID,stats45Win36);
% writetoFileTex('45Win60',pID,stats45Win60);
% writetoFileTex('60Win36',pID,stats60Win36);
% writetoFileTex('60Win60',pID,stats60Win60);
% writetoFileTex('120Win60',pID,stats120Win60);
%%write to csv python
% 
% writetoFilePy('120minWindow_avg_3stdev',pID,stats30Win24,30)
% writetoFilePy('120minWindow_avg_3stdev',pID,stats45Win24,45)
% 
% writetoFilePy('180minWindow_avg_3stdev',pID,stats30Win36,30)
% writetoFilePy('180minWindow_avg_3stdev',pID,stats45Win36,45)
% writetoFilePy('180minWindow_avg_3stdev',pID,stats60Win36,60)
% 
% writetoFilePy('300minWindow_avg_3stdev',pID,stats30Win60,30)
% writetoFilePy('300minWindow_avg_3stdev',pID,stats45Win60,45)
% writetoFilePy('300minWindow_avg_3stdev',pID,stats60Win60,60)
% writetoFilePy('300minWindow_avg_3stdev',pID,stats120Win60,120)

% % %% write to txt for python
% writetoFilePy('delta30Win24',pID,stats30Win24,30,120)
% writetoFilePy('delta45Win24',pID,stats45Win24,45,120)
% 
% writetoFilePy('delta30Win36',pID,stats30Win36,30,180)
% writetoFilePy('delta45Win36',pID,stats45Win36,45,180)
% writetoFilePy('delta60Win36',pID,stats60Win36,60,180)
% 
% writetoFilePy('delta30Win60',pID,stats30Win60,30,300)
% writetoFilePy('delta45Win60',pID,stats45Win60,45,300)
% writetoFilePy('delta60Win60',pID,stats60Win60,60,300)
% writetoFilePy('delta120Win60',pID,stats120Win60,120,300)
% 
% writetoFilePy('120minWindow',pID,stats30Win24,30,120)
% writetoFilePy('120minWindow',pID,stats45Win24,45,120)
% 
% writetoFilePy('180minWindow',pID,stats30Win36,30,180)
% writetoFilePy('180minWindow',pID,stats45Win36,45,180)
% writetoFilePy('180minWindow',pID,stats60Win36,60,180)
% 
% writetoFilePy('300minWindow',pID,stats30Win60,30,300)
% writetoFilePy('300minWindow',pID,stats45Win60,45,300)
% writetoFilePy('300minWindow',pID,stats60Win60,60,300)
% writetoFilePy('300minWindow',pID,stats120Win60,120,300)
%%

% writetoFilePyAll('120minWindow_3stdev',pID,modelFits30Win24,stats30Win24,30)
% writetoFilePyAll('120minWindow_3stdev',pID,modelFits45Win24,stats45Win24,45)
% 
% writetoFilePyAll('180minWindow_3stdev',pID,modelFits30Win36,stats30Win36,30)
% writetoFilePyAll('180minWindow_3stdev',pID,modelFits45Win36,stats45Win36,45)
% writetoFilePyAll('180minWindow_3stdev',pID,modelFits60Win36,stats60Win36,60)
% 
% writetoFilePyAll('300minWindow_3stdev',pID,modelFits30Win60,stats30Win60,30)
% writetoFilePyAll('300minWindow_3stdev',pID,modelFits45Win60,stats45Win60,45)
% writetoFilePyAll('300minWindow_3stdev',pID,modelFits60Win60,stats60Win60,60)
% writetoFilePyAll('300minWindow_3stdev',pID,modelFits120Win60,stats120Win60,120)

%%
close all
clear corrData corrStats;

corrData(MAX).delta6CGM=[];
corrData(MAX).delta6IOB=[];
corrData(MAX).delta9CGM=[];
corrData(MAX).delta12CGM=[];
corrData(MAX).delta24CGM=[];

corrStats(MAX).lags=[];
corrStats(MAX).bounds=[];
corrStats(MAX).order=[];
corrStats(MAX).value=[];
corrStats(MAX).top10=[];

RANGE=1:MAX;
RANGE(7)=[];
%take parcorr (partial autocorrelation) to determine how many terms are
%needed
for i=RANGE
    points6=1:6:length(patient(i).gtimes);
    corrData(i).delta6CGM=patient(i).gCGM(points6);
    corrData(i).delta6IOB=patient(i).gIOB(points6);
    points9=1:9:length(patient(i).gtimes);
    points12=1:12:length(patient(i).gtimes);
    points24=1:24:length(patient(i).gtimes);
    corrData(i).delta9CGM=patient(i).gCGM(points9);
    corrData(i).delta12CGM=patient(i).gCGM(points12);
    corrData(i).delta24CGM=patient(i).gCGM(points24);
end

for i=RANGE
    clear B I
    [corrStats(i).lags, ~, corrStats(i).bounds]=parcorr(patient(i).gCGM,min(length(patient(i).gCGM)-1,60));
    [B, I]= sort(corrStats(i).lags,'descend');
    corrStats(i).order=I;
    corrStats(i).value=B;
    corrStats(i).top10=I(1:min(10,length(corrStats(i).lags)));
    
    clear B I
    [AutocorrStats(i).lags, ~, AutocorrStats(i).bounds]=autocorr(patient(i).gCGM,min(length(patient(i).gCGM)-1,60));
    [B, I]= sort(AutocorrStats(i).lags,'descend');
    AutocorrStats(i).order=I;
    AutocorrStats(i).value=B;
    AutocorrStats(i).top10=I(1:min(10,length(AutocorrStats(i).lags)));
    
    clear B I
    [corrStats(i).lags6, ~, corrStats(i).bounds6]=parcorr(corrData(i).delta6CGM,min(length(corrData(i).delta6CGM)-1,10));
    [B, I]= sort(corrStats(i).lags6,'descend');
    corrStats(i).order_6=I;
    corrStats(i).value_6=B;
    corrStats(i).top10_6=I(1:min(10,length(corrStats(i).lags6)));
    
    clear B I
    [corrStats(i).lags9, ~, corrStats(i).bounds9]=parcorr(corrData(i).delta9CGM,min(length(corrData(i).delta9CGM)-1,7));
    [B, I]= sort(corrStats(i).lags9,'descend');
    corrStats(i).order_9=I;
    corrStats(i).value_9=B;
    corrStats(i).top10_9=I(1:min(10,length(corrStats(i).lags9)));
    
    clear B I
    [corrStats(i).lags12, ~, corrStats(i).bounds12]=parcorr(corrData(i).delta12CGM,min(length(corrData(i).delta12CGM)-1,5));
    [B, I]= sort(corrStats(i).lags12,'descend');
    corrStats(i).order_12=I;
    corrStats(i).value_12=B;
    corrStats(i).top10_12=I(1:min(10,length(corrStats(i).lags12)));

    clear B I
    [corrStats(i).lags24, ~, corrStats(i).bounds24]=parcorr(corrData(i).delta24CGM,min(length(corrData(i).delta24CGM)-1,3));
    [B, I]= sort(corrStats(i).lags24,'descend');
    corrStats(i).order_24=I;
    corrStats(i).value_24=B;
    corrStats(i).top10_24=I(1:min(10,length(corrStats(i).lags24)));
    
    figure(i)
    subplot(5,1,1)
    parcorr(corrData(i).delta6CGM,min(length(corrData(i).delta6CGM)-1,10))
    xlabel('CGM delta6')
    subplot(5,1,2)
    parcorr(corrData(i).delta9CGM,min(length(corrData(i).delta9CGM)-1,7))
	xlabel('CGM delta9')
    subplot(5,1,3)
    parcorr(corrData(i).delta12CGM,min(length(corrData(i).delta12CGM)-1,5))
	xlabel('CGM delta12')
    subplot(5,1,4)
    parcorr(corrData(i).delta24CGM,min(length(corrData(i).delta24CGM)-1,3))
	xlabel('CGM delta24')
    subplot(5,1,5)
    parcorr(patient(i).gCGM,min(length(patient(i).gCGM)-1,60))
    xlabel('CGM delta0')
    
    figure(100+i)
    autocorr(patient(i).gCGM,min(length(patient(i).gCGM)-1))
    xlabel('CGM delta0')
    ylabel('autocorrelation')
    title(strcat({'PSO3-001-0001 session ID: '}, num2str(cleanpatient(i).SessionID(1))))
end
%%
topFits=[corrStats.top10];
MEANS.all=mean(topFits,2);
[MODES.all(:,1),MODES.all(:,2)]=mode(topFits,2);
AA=horzcat(corrStats.order);
BB=horzcat(corrStats.value);
MODES.all(1,3)=mean(BB(1,top1));
top2=find(AA(2,:)==2);
MODES.all(1,4)=mean(BB(2,top2));

topFits6=[corrStats.top10_6];
MEANS.w6=mean(topFits6,2);
[MODES.w6(:,1),MODES.w6(:,2)]=mode(topFits6,2);
AA=horzcat(corrStats.order_6);
top1=find(AA(1,:)==1);
BB=horzcat(corrStats.value_6);
MODES.w6(1,3)=mean(BB(1,top1));
top2=find(AA(2,:)==2);
MODES.w6(1,4)=mean(BB(2,top2));

topFits9=[corrStats.top10_9];
MEANS.w9=mean(topFits9,2);
[MODES.w9(:,1),MODES.w9(:,2)]=mode(topFits9,2);
AA=horzcat(corrStats.order_9);
top1=find(AA(1,:)==1);
BB=horzcat(corrStats.value_9);
MODES.w9(1,3)=mean(BB(1,top1));
top2=find(AA(2,:)==2);
MODES.w9(1,4)=mean(BB(2,top2));

topFits12=[corrStats.top10_12];
MEANS.w12=mean(topFits12,2);
[MODES.w12(:,1),MODES.w12(:,2)]=mode(topFits12,2);
AA=horzcat(corrStats.order_12);
top1=find(AA(1,:)==1);
BB=horzcat(corrStats.value_12);
MODES.w12(1,3)=mean(BB(1,top1));
top2=find(AA(2,:)==2);
MODES.w12(1,4)=mean(BB(2,top2));

AtopFits=[AutocorrStats.top10];
MEANS.auto=mean(AtopFits,2);
[MODES.auto(:,1),MODES.auto(:,2)]=mode(AtopFits,2);

%take autocorr to determine time lags
% for i=1:MAX2
%     figure(100+i)
%     subplot(2,1,1)
%     autocorr(patient(i).gCGM,25)
%     xlabel('CGM')
%     subplot(2,1,2)
%     autocorr(patient(i).gIOB,25)
%     xlabel('IOB')
% end

%%
PeakCorr(length(RANGE)).pks=[];
PeakCorr(length(RANGE)).locs=[];
PeakCorr(length(RANGE)).lags=[];

Peaks(length(RANGE)).CGMpks=[];
n=4;
for i=RANGE
    lagged=(patient(i).gCGM(n+1:end)-patient(i).gCGM(1:end-n))/5;
    %smooth patient data
    smoothed=conv(patient(i).gCGM,[0.25, .5, 0.25],'valid');

    [acor, lag]=xcorr(-1*patient(i).gIOB(n+1:end), smooth(lagged));
    figure(300+i)
    subplot(4,1,1)
    plot(patient(i).gtimes(2:end-1),patient(i).gCGM(2:end-1),patient(i).gtimes(2:end-1),smoothed)
    xlabel('time since start trial')
    ylabel('CGM')
    subplot(4,1,2)
    plot(patient(i).gtimes,patient(i).gIOB)
    xlabel('time since start trial')
    ylabel('IOB')
    subplot(4,1,3)
    smoothX=smooth(lag,acor,0.1,'rloess');
    plot(lag,smoothX)
    xlabel('Lag')
    ylabel('smoothed cross correlation')
    subplot(4,1,4)
    plot(patient(i).gtimes(n+1:end),lagged,patient(i).gtimes(n+1:end),smooth(lagged))
    
    [PeakCorr(i).pks,PeakCorr(i).locs]=findpeaks(smoothX);
    PeakCorr(i).lags=lag(PeakCorr(i).locs);
    PeakCorr(i).diffLags=diff(PeakCorr(i).lags);
    PeakCorr(i).minLags=min(abs(PeakCorr(i).lags));
    
    [Peaks(i).CGMpks,Peaks(i).CGMlocs]=findpeaks(patient(i).gCGM);
    [Peaks(i).IOBpks,Peaks(i).IOBlocs]=findpeaks(patient(i).gIOB);
end

%% cluster the peaks
allpeaks=[PeakCorr.minLags]';
opts = statset('Display','final');
[idx,C] = kmeans(allpeaks,2,'Distance','cityblock',...
    'Replicates',5,'Options',opts);
%%
figure(400);
histogram(allpeaks(idx==1),40)
hold on
histogram(allpeaks(idx==2),40)
plot(C,[sum(idx==1)/10,sum(idx==2)/10],'kx',...
     'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids for Min Peaks'
hold off