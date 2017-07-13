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
    modelFits30min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits30min(i).predict(m)=patient(i).gCGM(m);
    end
    for m = delta:(tEND-gDelta)
        modelFits30min(i).predict(m+gDelta)=modelFits30min(i).mean(1)*patient(i).gCGM(m)+modelFits30min(i).mean(2)*patient(i).gCGM(m-gDelta)+modelFits30min(i).mean(3)*patient(i).gIOB(m-iDelta);
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

gDelta=24;
iDelta=24;
delta=max(gDelta,iDelta)+1;
for i=ID
    tEND=length(patient(i).gtimes);
    modelFits120min(i).predict=zeros(length(patient(i).gCGM),1);
    if tEND>(delta+gDelta+1)
        
        for m=1:delta+gDelta
            modelFits120min(i).predict(m)=patient(i).gCGM(m);
        end
        for m=delta:tEND-gDelta
            modelFits120min(i).predict(m+gDelta)=modelFits120min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits120min(i).mean(2)*patient(i).gCGM(m)+modelFits120min(i).mean(3)*patient(i).gIOB(m-iDelta);
        end
    else
        for m=1:tEND
            modelFits120min(i).predict(m)=patient(i).gCGM(m);
        end
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

%%
for i=1:MAX
    figure(200+i)
    subplot(2,2,1)
    hist(modelFits30min(i).predict(13:end)-patient(i).gCGM(13:end),200)
    title('Residuals for averaged parameters 30min')
    subplot(2,2,2)
    hist(modelFits45min(i).predict(19:end)-patient(i).gCGM(19:end),200)
    title('Residuals for averaged parameters 45min')    
    subplot(2,2,3)
    hist(modelFits60min(i).predict(25:end)-patient(i).gCGM(25:end),200)
    title('Residuals for averaged parameters 60min')    
    subplot(2,2,4)
    hist(modelFits120min(i).predict(48:end)-patient(i).gCGM(48:end),200)
    title('Residuals for averaged parameters 120min')
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
    for m = delta:tEND-gDelta
        modelFits30min(i).predict(m+gDelta)=modelFits30min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits30min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits30min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
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
    for m = delta:tEND-gDelta
        modelFits45min(i).predict(m+gDelta)=modelFits45min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits45min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits45min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
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
    for m = delta:tEND-gDelta
        modelFits60min(i).predict(m+gDelta)=modelFits60min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits60min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits60min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
    end
end

gDelta=24;
iDelta=24;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    modelFits120min(i).predict=zeros(length(patient(i).gCGM),1);
    if tEND>(delta+gDelta+1)
        for m=1:delta+gDelta+1
            modelFits120min(i).predict(m)=patient(i).gCGM(m);
        end
        for m = delta:tEND-gDelta
            modelFits120min(i).predict(m+gDelta)=modelFits120min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits120min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits120min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
        end
    else
        for m=1:tEND
            modelFits120min(i).predict(m)=patient(i).gCGM(m);
        end
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

%%

%% plot fit models 30min

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
    for m = delta:tEND
        modelFits30min(i).predict(m)=modelFits30min(i).mean(1)*patient(i).gCGM(m)+modelFits30min(i).mean(2)*patient(i).gCGM(m-gDelta)+modelFits30min(i).mean(3)*patient(i).gIOB(m-iDelta);
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
    for m=delta:tEND
        modelFits60min(i).predict(m)=modelFits60min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits60min(i).mean(2)*patient(i).gCGM(m)+modelFits60min(i).mean(3)*patient(i).gIOB(m-iDelta);
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
    for m=delta:tEND
        modelFits45min(i).predict(m)=modelFits45min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits45min(i).mean(2)*patient(i).gCGM(m)+modelFits45min(i).mean(3)*patient(i).gIOB(m-iDelta);
    end
end

gDelta=24;
iDelta=24;
delta=max(gDelta,iDelta)+1;
for i=ID
    tEND=length(patient(i).gtimes);
    modelFits120min(i).predict=zeros(length(patient(i).gCGM),1);
    if tEND>(delta+gDelta+1)
        
        for m=1:delta+gDelta
            modelFits120min(i).predict(m)=patient(i).gCGM(m);
        end
        for m=delta:tEND
            modelFits120min(i).predict(m)=modelFits120min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits120min(i).mean(2)*patient(i).gCGM(m)+modelFits120min(i).mean(3)*patient(i).gIOB(m-iDelta);
        end
    else
        for m=1:tEND
            modelFits120min(i).predict(m)=patient(i).gCGM(m);
        end
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

%%

for i=1:5
figure(100+i)
subplot(2,2,1)
hist(modelFits30min(i).RES,200)
title('Residuals 30min')
subplot(2,2,2)
hist(modelFits45min(i).RES,200)
title('Residuals 45min')
subplot(2,2,3)
hist(modelFits60min(i).RES,200)
title('Residuals 60min')
subplot(2,2,4)
hist(modelFits120min(i).RES,200)
title('Residuals 120min')
end

%% GLOBAL OPTIMIZATION
MAX=10;

a0=[.5 .5 -2.5];
lb = [-2 -2 -10];
ub=[4 4 0];
gDelta=6;
iDelta=6;
tic
[modelFits30min, stats30min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc
%%
a0=[.5 .5 -2.5];
lb = [-2 -2 -10];
ub=[4 4 0];
gDelta=6;
iDelta=6;
tic
[modelFits30minGlob, stats30minGlob]=GlobalRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc

%% Plots
ID=1:20;
gDelta=6;
iDelta=6;
delta=max(gDelta,iDelta)+1;


for i=1:MAX
    tEND=length(patient(i).gtimes);
    modelFits30min(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits30min(i).predict(m)=patient(i).gCGM(m);
    end
    for m = delta:tEND-gDelta
        modelFits30min(i).predict(m+gDelta)=stats30min.mean(1)*patient(i).gCGM(m-gDelta)+stats30min.mean(2)*patient(i).gCGM(m)+stats30min.mean(3)*patient(i).gIOB(m-iDelta);
    end
end

for i=1:MAX
    tEND=length(patient(i).gtimes);
    modelFits30minGlob(i).predict=zeros(length(patient(i).gCGM),1);
    for m=1:delta+gDelta+1
        modelFits30minGlob(i).predict(m)=patient(i).gCGM(m);
    end
    for m = delta:tEND-gDelta
        modelFits30minGlob(i).predict(m+gDelta)=stats30minGlob.mean(1)*patient(i).gCGM(m-gDelta)+stats30minGlob.mean(2)*patient(i).gCGM(m)+stats30minGlob.mean(3)*patient(i).gIOB(m-iDelta);
    end
end

for i=1:MAX
    figure(i)
    subplot(2,1,1)
    plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,modelFits30min(i).predict,patient(i).gtimes,modelFits30minGlob(i).predict)
    legend('CGM data','30min','30minGlob')
    title('CGM actual vs Predicted values for various time steps')
    ylabel('Glucose Measurements (mg/dL)')
    xlabel('minutes since start of trial')
    subplot(2,1,2)
    plot(patient(i).times,patient(i).IOB)
    ylabel('Insulin On Board')
    xlabel('minutes since start of trial')
    title(allfiles(i).name)
end

for i=1:MAX
figure(100+i)
subplot(2,1,1)
hist(modelFits30min(i).predict(13:end)-patient(i).gCGM(13:end),200)
title('Residuals 30min')
subplot(2,1,2)
hist(modelFits30minGlob(i).predict(13:end)-patient(i).gCGM(13:end),200)
title('Residuals 30min Global')
end
%%
for i=1:MAX
    avg30(i).res=modelFits30min(i).predict(13:end)-patient(i).gCGM(13:end);
    glob30(i).res=modelFits30minGlob(i).predict(13:end)-patient(i).gCGM(13:end);
end

avg30mean=padcat(avg30(1:end).res);
glob30mean=padcat(glob30(1:end).res);
statsAvg30.mean=nanmean(nanmean(abs(avg30mean)));
statsGlob30.mean=nanmean(nanmean(abs(glob30mean)));
statsAvg30.stdev=nanstd(nanstd(avg30mean));
statsGlob30.stdev=nanstd(nanstd(glob30mean));
statsAvg30.min=min(min(abs(avg30mean)));
statsGlob30.min=min(min(abs(glob30mean)));
statsAvg30.max=max(max(abs(avg30mean)));
statsGlob30.max=max(max(abs(glob30mean)));
statsAvg30.c95=1.96*statsAvg30.stdev;
statsGlob30.c95=1.96*statsGlob30.stdev;


