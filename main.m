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
MAX2=numel(cleanpatient);
a0=[.6 .4 -2.5];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=6;
Delta1=6;
stepsz=24;
ovlp=5;

[modelFits30Win24, stats30Win24] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);
%[modelFits30Win24_dirty, stats30Win24_dirty]=WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

stepsz=36;
ovlp=5;
[modelFits30Win36, stats30Win36] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);

stepsz=60;
ovlp=5;
[modelFits30Win60, stats30Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);

% 45min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=9;
Delta1=9;

stepsz=24;
ovlp=5;
[modelFits45Win24, stats45Win24] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);


stepsz=36;
ovlp=5;
[modelFits45Win36, stats45Win36] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits45Win60, stats45Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);

% 60min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=12;
Delta1=12;

stepsz=36;
ovlp=5;
[modelFits60Win36, stats60Win36] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits60Win60, stats60Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);

%120min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
Delta1=24;
Delta1=24;

stepsz=60;
ovlp=5;

[modelFits120Win60, stats120Win60] = WindowRegModelFit(a0,lb,ub,Delta1,Delta1,cleanpatient,MAX2,stepsz,ovlp);
toc
%% plots for presentation
close all
Delta1=6;
Delta2=9;

for i=2:12
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
        predictMAX(n)=max(predictmax1(n),predictmax2(n-diffdelta));
        predictMIN(n)=min(predictmin1(n),predictmin2(n-diffdelta));
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
    title(strcat('PSO3-001-0001 trial',num2str(i)),'FontSize',16)
    subplot(2,1,2)
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gIOB(2*Delta1+1:end))
    xlabel('time since start trial','FontSize',14)
    ylabel('IOB (units)','FontSize',14)
end
%%
close all
for i=2:15
    A1=[patient(i).gCGM(1:end-2*Delta1), patient(i).gCGM(Delta1+1:end-Delta1), patient(i).gIOB(1:end-2*Delta1)];
    A2=[patient(i).gCGM(1:end-2*Delta2), patient(i).gCGM(Delta2+1:end-Delta2), patient(i).gIOB(1:end-2*Delta2)];
    
    %x=[modelFits30Win36(i).MEAN(1); modelFits30Win36(i).MEAN(2); modelFits30Win36(i).MEAN(3)];
    x1=[stats30Win36.mean(1); stats30Win36.mean(2); stats30Win36.mean(3)];
    x2=[stats45Win36.mean(1); stats45Win36.mean(2); stats45Win36.mean(3)];
    predict1=A1*x1;
    predict2=A2*x2;
    predictmin1=predict1-stats30Win36.RES95;
    predictmax1=predict1+stats30Win36.RES95;
    
    predictmin2=predict2-stats45Win36.RES95;
    predictmax2=predict2+stats45Win36.RES95;
    
    predictMAX=zeros(size(predict1));
    predictMIN=zeros(size(predict1));
    diffdelta=2*Delta2-2*Delta1;
    for n=1:diffdelta
        predictMAX(n)=predictmax1(n);
        predictMIN(n)=predictmin1(n);
    end
    
    for n=diffdelta+1:length(predictMAX)
        predictMAX(n)=max(predictmax1(n),predictmax2(n-diffdelta));
        predictMIN(n)=min(predictmin1(n),predictmin2(n-diffdelta));
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
    title(strcat('PSO3-001-0001 trial',num2str(i)),'FontSize',16)
    subplot(2,1,2)
    plot(patient(i).gtimes(2*Delta1+1:end),patient(i).gIOB(2*Delta1+1:end))
    xlabel('time since start trial','FontSize',14)
    ylabel('IOB (units)','FontSize',14)
end
%%
close all
for i=2:12
    A1=[patient(i).gCGM(1:end-2*Delta1), patient(i).gCGM(Delta1+1:end-Delta1), patient(i).gIOB(1:end-2*Delta1)];
    A2=[patient(i).gCGM(1:end-2*Delta2), patient(i).gCGM(Delta2+1:end-Delta2), patient(i).gIOB(1:end-2*Delta2)];
    
    %x=[modelFits30Win60(i).MEAN(1); modelFits30Win60(i).MEAN(2); modelFits30Win60(i).MEAN(3)];
    x1=[stats30Win60.mean(1); stats30Win60.mean(2); stats30Win60.mean(3)];
    x2=[stats45Win60.mean(1); stats45Win60.mean(2); stats45Win60.mean(3)];
    predict1=A1*x1;
    predict2=A2*x2;
    predictmin1=predict1-stats30Win60.RES95;
    predictmax1=predict1+stats30Win60.RES95;
    
    predictmin2=predict2-stats45Win60.RES95;
    predictmax2=predict2+stats45Win60.RES95;
    
    predictMAX=zeros(size(predict1));
    predictMIN=zeros(size(predict1));
    diffdelta=2*Delta2-2*Delta1;
    for n=1:diffdelta
        predictMAX(n)=predictmax1(n);
        predictMIN(n)=predictmin1(n);
    end
    
    for n=diffdelta+1:length(predictMAX)
        predictMAX(n)=max(predictmax1(n),predictmax2(n-diffdelta));
        predictMIN(n)=min(predictmin1(n),predictmin2(n-diffdelta));
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
    title(strcat('PSO3-001-0001 trial',num2str(i)),'FontSize',16)
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

