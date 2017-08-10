%30 June 2017
%Taisa Kushner
%fits for data

%turn off the beeping
beep off 
clear all

%pull all files for given patient
pID='PSO3-001-0001';

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

%% Optimization by window
tic 
a0=[.6 .4 -2.5];
lb = [0 0 -100];
ub=[1.5 1.5 0];
gDelta=6;
iDelta=6;
stepsz=24;
ovlp=5;

[modelFits30Win24, stats30Win24] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=36;
ovlp=5;
[modelFits30Win36, stats30Win36] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

stepsz=60;
ovlp=5;
[modelFits30Win60, stats30Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

% 45min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
gDelta=9;
iDelta=9;

stepsz=24;
ovlp=5;
[modelFits45Win24, stats45Win24] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=36;
ovlp=5;
[modelFits45Win36, stats45Win36] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits45Win60, stats45Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

% 60min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
gDelta=12;
iDelta=12;

stepsz=36;
ovlp=5;
[modelFits60Win36, stats60Win36] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits60Win60, stats60Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

%120min

a0=[.6 .4 -10];
lb = [0 0 -100];
ub=[1.5 1.5 0];
gDelta=24;
iDelta=24;

stepsz=60;
ovlp=5;

[modelFits120Win60, stats120Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);
toc

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
% %% write to txt for python
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

writetoFilePy('120minWindow',pID,stats30Win24,30,120)
writetoFilePy('120minWindow',pID,stats45Win24,45,120)

writetoFilePy('180minWindow',pID,stats30Win36,30,180)
writetoFilePy('180minWindow',pID,stats45Win36,45,180)
writetoFilePy('180minWindow',pID,stats60Win36,60,180)

writetoFilePy('300minWindow',pID,stats30Win60,30,300)
writetoFilePy('300minWindow',pID,stats45Win60,45,300)
writetoFilePy('300minWindow',pID,stats60Win60,60,300)
writetoFilePy('300minWindow',pID,stats120Win60,120,300)
%%

% writetoFilePyAll('120minWindow',pID,modelFits30Win24,stats30Win24,30)
% writetoFilePyAll('120minWindow',pID,modelFits45Win24,stats45Win24,45)
% 
% writetoFilePyAll('180minWindow',pID,modelFits30Win36,stats30Win36,30)
% writetoFilePyAll('180minWindow',pID,modelFits45Win36,stats45Win36,45)
% writetoFilePyAll('180minWindow',pID,modelFits60Win36,stats60Win36,60)
% 
% writetoFilePyAll('300minWindow',pID,modelFits30Win60,stats30Win60,30)
% writetoFilePyAll('300minWindow',pID,modelFits45Win60,stats45Win60,45)
% writetoFilePyAll('300minWindow',pID,modelFits60Win60,stats60Win60,60)
% writetoFilePyAll('300minWindow',pID,modelFits120Win60,stats120Win60,120)

