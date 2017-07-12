%30 June 2017
%compiled fits to pull function instead of all together

%pull all files for given patient
allfiles=dir('../outputs/byPatient/session-PSO3-001-*');
%how many are there
numfiles=size(allfiles);
MAX=20;

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

tic
parfor i=1:MAX
    [patient(i).SessionID,patient(i).datetime,patient(i).CGM,patient(i).IOB,patient(i).Bolus,patient(i).BkgInsulin] = ...
        importCGMDATA(strcat('../outputs/byPatient/',allfiles(i).name));
    
    formatin='HH:MM:SS mm/dd/yyyy';
    patient(i).times=datenum(patient(i).datetime,formatin);
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
    
end
toc

%% minimize least squares for 30min
a0=[.8 .2 -.4]; 
lb = [0 0 -5];
ub=[4 4 0];
gDelta=6;
iDelta=6;
tic
[modelFits30min, stats30min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc

%% minimize least squares for 45min
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[4 4 0];
gDelta=9;
iDelta=9;
tic
[modelFits45min, stats45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc
%% minimize least squares for 60min
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[4 4 0];
gDelta=12;
iDelta=12;
tic
[modelFits60min, stats60min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc
%% 120min
a0=[2 1 -10]; 
lb = [0 0 -50];
ub=[4 4 0];
gDelta=24;
iDelta=24;
tic
[modelFits120min, stats120min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc

%% mixed 120 & 45
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[4 4 0];
gDelta=24;
iDelta=9;
tic
[modelFitsG120I45min, statsG120I45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc
%mixed 60 & 45
a0=[2 1 -10]; 
lb = [0 0 -30];
ub=[4 4 0];
gDelta=12;
iDelta=9;
tic
[modelFitsG60I45min, statsG60I45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc

%% plot fit models 30min

ID=1:20;
gDelta=6;
iDelta=6;
delta=max(gDelta,iDelta)+1;


for i=ID
    tEND=length(patient(i).gtimes);
    modelFits30min(i).predictAVG=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits30min(i).predictAVG(m)=patient(i).gCGM(m);
    end
    for m = delta:(tEND-gDelta)
        modelFits30min(i).predictAVG(m+gDelta)=modelFits30min(i).mean(1)*patient(i).gCGM(m)+modelFits30min(i).mean(2)*patient(i).gCGM(m-gDelta)+modelFits30min(i).mean(3)*patient(i).gIOB(m-iDelta);
    end
end

gDelta=12;
iDelta=12;
delta=max(gDelta,iDelta)+1;
for i=ID
        tEND=length(patient(i).gtimes);
    modelFits60min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta
        modelFits60min(i).predict(m)=patient(i).gCGM(m);
    end
    for m=delta:tEND-gDelta
        modelFits60min(i).predict(m+gDelta)=modelFits60min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits60min(i).mean(2)*patient(i).gCGM(m)+modelFits60min(i).mean(3)*patient(i).gIOB(m-iDelta);
    end
end

gDelta=9;
iDelta=9;
delta=max(gDelta,iDelta)+1;
for i=ID
        tEND=length(patient(i).gtimes);
        modelFits45min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta
        modelFits45min(i).predict(m)=patient(i).gCGM(m);
    end
    for m=delta:tEND-gDelta
        modelFits45min(i).predict(m+gDelta)=modelFits45min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits45min(i).mean(2)*patient(i).gCGM(m)+modelFits45min(i).mean(3)*patient(i).gIOB(m-iDelta);
    end
end

for i=ID
    % plot
    figure(i)
    subplot(2,1,1)
    plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,modelFits30min(i).predictAVG,patient(i).gtimes,modelFits60min(i).predict,patient(i).gtimes,modelFits45min(i).predict)
    %plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,modelFits30min(i).predict,patient(i).gtimes,modelFits30min(i).predictAVG)
    ylabel('CGM')
    xlabel('minutes since start')
    title(strcat(patient(i).datetime(1)))
    subplot(2,1,2)
    plot(patient(i).times,patient(i).IOB)
    ylabel('IOB')
    xlabel('minutes since start')
    title(strcat(patient(i).datetime(end)))
end

%%

ID=1:20;
gDelta=6;
iDelta=6;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    modelFits30min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits30min(i).predict(m)=patient(i).gCGM(m);
    end
    n=1;
    for m = (delta+gDelta):tEND-3
        modelFits30min(i).predict(m)=modelFits30min(i).Fits(n,1)*patient(i).gCGM(n)+modelFits30min(i).Fits(n,2)*patient(i).gCGM(n+gDelta)+modelFits30min(i).Fits(n,3)*patient(i).gIOB(n+iDelta);
        n=n+1;
    end
    for m= tEND-3:tEND
        modelFits30min(i).predict(m)=patient(i).gCGM(m);
    end
end

gDelta=9;
iDelta=9;
delta=max(gDelta,iDelta)+1;


for i=ID
    tEND=length(patient(i).gtimes);
    modelFits45min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits45min(i).predict(m)=patient(i).gCGM(m);
    end
    for m = delta:tEND-3
        modelFits45min(i).predict(m+gDelta)=modelFits45min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits45min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits45min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
    end
    for m= tEND-3:tEND
        modelFits45min(i).predict(m)=patient(i).gCGM(m);
    end
end

gDelta=12;
iDelta=12;
delta=max(gDelta,iDelta)+1;


for i=ID
    tEND=length(patient(i).gtimes);
    modelFits60min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits60min(i).predict(m)=patient(i).gCGM(m);
    end
    for m = delta:tEND-3
        modelFits60min(i).predict(m+gDelta)=modelFits60min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits60min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits60min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
    end
    for m= tEND-3:tEND
        modelFits60min(i).predict(m)=patient(i).gCGM(m);
    end
end

gDelta=24;
iDelta=24;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    modelFits120min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits120min(i).predict(m)=patient(i).gCGM(m);
    end
    for m = delta:tEND-3
        modelFits120min(i).predict(m+gDelta)=modelFits120min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits120min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits120min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
    end
    for m= tEND-3:tEND
        modelFits120min(i).predict(m)=patient(i).gCGM(m);
    end
end

for i=ID
    % plot
    figure(i)
    subplot(2,1,1)
    plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,modelFits30min(i).predict,patient(i).gtimes,modelFits45min(i).predict,patient(i).gtimes,modelFits60min(i).predict,patient(i).gtimes,modelFits120min(i).predict)
    legend('CGM data','30min','45min','60min','120min')
    title('CGM actual vs Predicted values for various time steps')
    ylabel('Glucose Measurements (mg/dL)')
    xlabel('minutes since start of trial')
    subplot(2,1,2)
    plot(patient(i).times,patient(i).IOB)
    ylabel('Insulin On Board')
    xlabel('minutes since start of trial')
    title(allfiles(i).name)
end