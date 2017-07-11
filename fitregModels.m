%30 June 2017
%compiled fits to pull function instead of all together

%pull all files for given patient
allfiles=dir('./outputs/byPatient/session-PSO3-001-*');
%how many are there
numfiles=size(allfiles);
MAX=657;
%%

%Initialize patient struct for speed
patient(MAX).SessionID=[];
patient(MAX).datetime=[];
patient(MAX).CGM=[];
patient(MAX).IOB=[];
patient(MAX).Bolus=[];
patient(MAX).BkgInsulin=[];
patient(MAX).times=[];
patient(MAX).gtimes=[];
patient(MAX).gCGM=[];
patient(MAX).gIOB=[];

for i=1:MAX
    [patient(i).SessionID,patient(i).datetime,patient(i).CGM,patient(i).IOB,patient(i).Bolus,patient(i).BkgInsulin] = ...
        importCGMDATA(strcat('./Outputs/byPatient/',allfiles(i).name));
    
    patient(i).times=datenum(patient(i).datetime);
    %scaling step=5min timescale - same for all
    step=patient(i).times(3)-patient(i).times(2);
    %rescale times based on start time
    patient(i).times=patient(i).times-patient(i).times(1);
    patient(i).times=patient(i).times./step*5;
    
    %pull only times after first 180min burn-in
    gindex=find((patient(i).times>=180)-ismissing(patient(i).IOB));
    patient(i).gtimes=patient(i).times(gindex);
    patient(i).gCGM=patient(i).CGM(gindex);
    patient(i).gIOB=patient(i).IOB(gindex);
    
end

%% minimize least squares for 30min
a0=[.8 .2 -.4]; 
lb = [0 0 -5];
ub=[1 2 0];
gDelta=6;
iDelta=6;
delta=max(gDelta,iDelta)+1;

[modelFits30min, stats30min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
%% minimize least squares for 45min
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[1 2 0];
gDelta=9;
iDelta=9;
delta=max(gDelta,iDelta)+1;

[modelFits45min, stats45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
%% minimize least squares for 60min
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[2 2 0];
gDelta=12;
iDelta=12;
delta=max(gDelta,iDelta)+1;

[modelFits60min, stats60min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);

%%120min
a0=[2 1 -10]; 
lb = [0 0 -50];
ub=[2 2 0];
gDelta=24;
iDelta=24;
delta=max(gDelta,iDelta)+1;
[modelFits120min, stats120min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);


%%mixed 120 & 45
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[2 2 0];
gDelta=24;
iDelta=9;
delta=max(gDelta,iDelta)+1;
[modelFitsG120I45min, statsG120I45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);

%mixed 60 & 45
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[2 2 0];
gDelta=12;
iDelta=9;
delta=max(gDelta,iDelta)+1;
[modelFitsG60I45min, statsG60I45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);