%30 June 2017
%compiled fits to pull function instead of all together
%13 July added global optimization scheme to bottom

%pull all files for given patient
allfiles=dir('../outputs/byPatient/session-PSO3-001-*');
%how many are there
numfiles=size(allfiles);
MAX=657;

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
    
    if isempty(patient(i).gtimes)
        patient(i).gtimes=[0,0];
        patient(i).gCGM=0;
        patient(i).gIOB=0;
    end
end
toc

%% minimize least squares for 30min
a0=[.8 .2 -.4];
lb = [0 0 -5];
ub=[2 2 0];
gDelta=6;
iDelta=6;
tic
[modelFits30min, stats30min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc

% %% minimize least squares for 45min
% a0=[2 1 -10];
% lb = [0 0 -30];
% ub=[4 4 0];
% gDelta=9;
% iDelta=9;
% tic
% [modelFits45min, stats45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
% toc
% %% minimize least squares for 60min
% a0=[2 1 -10];
% lb = [0 0 -30];
% ub=[4 4 0];
% gDelta=12;
% iDelta=12;
% tic
% [modelFits60min, stats60min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
% toc
% %% 120min
% a0=[2 1 -10];
% lb = [0 0 -50];
% ub=[4 4 0];
% gDelta=24;
% iDelta=24;
% tic
% [modelFits120min, stats120min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
% toc
% 
% %% mixed 120 & 45
% a0=[2 1 -10];
% lb = [0 0 -30];
% ub=[4 4 0];
% gDelta=24;
% iDelta=9;
% tic
% [modelFitsG120I45min, statsG120I45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
% toc
% %mixed 60 & 45
% a0=[2 1 -10];
% lb = [0 0 -30];
% ub=[4 4 0];
% gDelta=12;
% iDelta=9;
% tic
% [modelFitsG60I45min, statsG60I45min]=RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
% toc
% 
% %% plot fit models for predictions using averaged param by patient trial
% 
% ID=1:MAX;
% gDelta=6;
% iDelta=6;
% delta=max(gDelta,iDelta)+1;
% 
% 
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits30min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta+1
%         modelFits30min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m = delta:(tEND-gDelta)
%         modelFits30min(i).predict(m+gDelta)=modelFits30min(i).mean(1)*patient(i).gCGM(m)+modelFits30min(i).mean(2)*patient(i).gCGM(m-gDelta)+modelFits30min(i).mean(3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=12;
% iDelta=12;
% delta=max(gDelta,iDelta)+1;
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits60min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta
%         modelFits60min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m=delta:tEND-gDelta
%         modelFits60min(i).predict(m+gDelta)=modelFits60min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits60min(i).mean(2)*patient(i).gCGM(m)+modelFits60min(i).mean(3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=9;
% iDelta=9;
% delta=max(gDelta,iDelta)+1;
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits45min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta
%         modelFits45min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m=delta:tEND-gDelta
%         modelFits45min(i).predict(m+gDelta)=modelFits45min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits45min(i).mean(2)*patient(i).gCGM(m)+modelFits45min(i).mean(3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=24;
% iDelta=24;
% delta=max(gDelta,iDelta)+1;
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits120min(i).predict=zeros(length(patient(i).gCGM),1);
%     if tEND>(delta+gDelta+1)
%         
%         for m=1:delta+gDelta
%             modelFits120min(i).predict(m)=patient(i).gCGM(m);
%         end
%         for m=delta:tEND-gDelta
%             modelFits120min(i).predict(m+gDelta)=modelFits120min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits120min(i).mean(2)*patient(i).gCGM(m)+modelFits120min(i).mean(3)*patient(i).gIOB(m-iDelta);
%         end
%     else
%         for m=1:tEND
%             modelFits120min(i).predict(m)=patient(i).gCGM(m);
%         end
%     end
% end
% 
% for i=ID
%     % plot
%     figure(i)
%     subplot(2,1,1)
%     plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,modelFits30min(i).predict,patient(i).gtimes,modelFits45min(i).predict,patient(i).gtimes,modelFits60min(i).predict,patient(i).gtimes,modelFits120min(i).predict)
%     legend('CGM data','30min','45min','60min','120min')
%     title('CGM actual vs Predicted values for various time steps')
%     ylabel('Glucose Measurements (mg/dL)')
%     xlabel('minutes since start of trial')
%     subplot(2,1,2)
%     plot(patient(i).times,patient(i).IOB)
%     ylabel('Insulin On Board')
%     xlabel('minutes since start of trial')
%     title(allfiles(i).name)
% end
% 
% %% plot histogram of residuals for above
% for i=1:20
%     figure(200+i)
%     subplot(2,2,1)
%     hist(modelFits30min(i).predict(13:end-13)-patient(i).gCGM(13:end-13),200)
%     title('Residuals for averaged parameters 30min')
%     subplot(2,2,2)
%     hist(modelFits45min(i).predict(19:end-19)-patient(i).gCGM(19:end-19),200)
%     title('Residuals for averaged parameters 45min')    
%     subplot(2,2,3)
%     hist(modelFits60min(i).predict(25:end-25)-patient(i).gCGM(25:end-25),200)
%     title('Residuals for averaged parameters 60min')    
%     subplot(2,2,4)
%     hist(modelFits120min(i).predict(48:end-48)-patient(i).gCGM(48:end-48),200)
%     title('Residuals for averaged parameters 120min')
% end
% %% plot fits using the found optimal value at each time window
% 
% ID=1:MAX;
% gDelta=6;
% iDelta=6;
% delta=max(gDelta,iDelta)+1;
% 
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits30min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta+1
%         modelFits30min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m = delta:tEND-gDelta
%         modelFits30min(i).predict(m+gDelta)=modelFits30min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits30min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits30min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
%     end
%     
% end
% 
% gDelta=9;
% iDelta=9;
% delta=max(gDelta,iDelta)+1;
% 
% 
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits45min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta+1
%         modelFits45min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m = delta:tEND-gDelta
%         modelFits45min(i).predict(m+gDelta)=modelFits45min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits45min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits45min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=12;
% iDelta=12;
% delta=max(gDelta,iDelta)+1;
% 
% 
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits60min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta+1
%         modelFits60min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m = delta:tEND-gDelta
%         modelFits60min(i).predict(m+gDelta)=modelFits60min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits60min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits60min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=24;
% iDelta=24;
% delta=max(gDelta,iDelta)+1;
% 
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits120min(i).predict=zeros(length(patient(i).gCGM),1);
%     if tEND>(delta+gDelta+1)
%         for m=1:delta+gDelta+1
%             modelFits120min(i).predict(m)=patient(i).gCGM(m);
%         end
%         for m = delta:tEND-gDelta
%             modelFits120min(i).predict(m+gDelta)=modelFits120min(i).Fits(m-gDelta,1)*patient(i).gCGM(m-gDelta)+modelFits120min(i).Fits(m-gDelta,2)*patient(i).gCGM(m)+modelFits120min(i).Fits(m-gDelta,3)*patient(i).gIOB(m-iDelta);
%         end
%     else
%         for m=1:tEND
%             modelFits120min(i).predict(m)=patient(i).gCGM(m);
%         end
%     end
% end
% 
% for i=ID
%     % plot
%     figure(i)
%     subplot(2,1,1)
%     plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,modelFits30min(i).predict,patient(i).gtimes,modelFits45min(i).predict,patient(i).gtimes,modelFits60min(i).predict,patient(i).gtimes,modelFits120min(i).predict)
%     legend('CGM data','30min','45min','60min','120min')
%     title('CGM actual vs Predicted values for various time steps')
%     ylabel('Glucose Measurements (mg/dL)')
%     xlabel('minutes since start of trial')
%     subplot(2,1,2)
%     plot(patient(i).times,patient(i).IOB)
%     ylabel('Insulin On Board')
%     xlabel('minutes since start of trial')
%     title(allfiles(i).name)
% end
%%
% 
% 
% %% plot fit models but not actually predicting (looks better....v weird and frustrating)
% 
% ID=1:20;
% gDelta=6;
% iDelta=6;
% delta=max(gDelta,iDelta)+1;
% 
% 
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits30min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta+1
%         modelFits30min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m = delta:tEND
%         modelFits30min(i).predict(m)=modelFits30min(i).mean(1)*patient(i).gCGM(m)+modelFits30min(i).mean(2)*patient(i).gCGM(m-gDelta)+modelFits30min(i).mean(3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=12;
% iDelta=12;
% delta=max(gDelta,iDelta)+1;
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits60min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta
%         modelFits60min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m=delta:tEND
%         modelFits60min(i).predict(m)=modelFits60min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits60min(i).mean(2)*patient(i).gCGM(m)+modelFits60min(i).mean(3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=9;
% iDelta=9;
% delta=max(gDelta,iDelta)+1;
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits45min(i).predict=zeros(length(patient(i).gCGM),1);
%     for m=1:delta+gDelta
%         modelFits45min(i).predict(m)=patient(i).gCGM(m);
%     end
%     for m=delta:tEND
%         modelFits45min(i).predict(m)=modelFits45min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits45min(i).mean(2)*patient(i).gCGM(m)+modelFits45min(i).mean(3)*patient(i).gIOB(m-iDelta);
%     end
% end
% 
% gDelta=24;
% iDelta=24;
% delta=max(gDelta,iDelta)+1;
% for i=ID
%     tEND=length(patient(i).gtimes);
%     modelFits120min(i).predict=zeros(length(patient(i).gCGM),1);
%     if tEND>(delta+gDelta+1)
%         
%         for m=1:delta+gDelta
%             modelFits120min(i).predict(m)=patient(i).gCGM(m);
%         end
%         for m=delta:tEND
%             modelFits120min(i).predict(m)=modelFits120min(i).mean(1)*patient(i).gCGM(m-gDelta)+modelFits120min(i).mean(2)*patient(i).gCGM(m)+modelFits120min(i).mean(3)*patient(i).gIOB(m-iDelta);
%         end
%     else
%         for m=1:tEND
%             modelFits120min(i).predict(m)=patient(i).gCGM(m);
%         end
%     end
% end
% 
% for i=ID
%     % plot
%     figure(i)
%     subplot(2,1,1)
%     plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,modelFits30min(i).predict,patient(i).gtimes,modelFits45min(i).predict,patient(i).gtimes,modelFits60min(i).predict,patient(i).gtimes,modelFits120min(i).predict)
%     legend('CGM data','30min','45min','60min','120min')
%     title('CGM actual vs Predicted values for various time steps')
%     ylabel('Glucose Measurements (mg/dL)')
%     xlabel('minutes since start of trial')
%     subplot(2,1,2)
%     plot(patient(i).times,patient(i).IOB)
%     ylabel('Insulin On Board')
%     xlabel('minutes since start of trial')
%     title(allfiles(i).name)
% end
% 
% %% Residuals for value optimized at each window (great)
% 
% for i=1:5
% figure(100+i)
% subplot(2,2,1)
% hist(modelFits30min(i).RES,200)
% title('Residuals 30min')
% subplot(2,2,2)
% hist(modelFits45min(i).RES,200)
% title('Residuals 45min')
% subplot(2,2,3)
% hist(modelFits60min(i).RES,200)
% title('Residuals 60min')
% subplot(2,2,4)
% hist(modelFits120min(i).RES,200)
% title('Residuals 120min')
% end
% 
%% GLOBAL OPTIMIZATION this is so much faster omg
%% minimize least squares for 30min
a0=[.5 .5 -2.5];
lb = [0 0 -10];
ub=[4 4 0];
gDelta=6;
iDelta=6;
tic
[GlobmodelFits30min, Globstats30min]=GlobalRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc

ID=1:MAX;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    if GlobmodelFits30min(i).RES >0 && ~isnan(GlobmodelFits30min(i).RES)
        GlobmodelFits30min(i).predict=[patient(i).gCGM(1:(delta+gDelta-1))' (GlobmodelFits30min(i).Fits(1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            GlobmodelFits30min(i).Fits(2)*patient(i).gCGM(delta:(end-gDelta))+ ...
            GlobmodelFits30min(i).Fits(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))']';
        
        GlobmodelFits30min(i).indivRes=abs(GlobmodelFits30min(i).predict-patient(i).gCGM);
    else
        GlobmodelFits30min(i).indivRes=0;
    end
end

Globstats30min.indivRes=vertcat(GlobmodelFits30min(1:end).indivRes);
Globstats30min.irMean=nanmean(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.irStdev=nanstd(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.irMin=min(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.irMax=max(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.ir95=1.96*Globstats30min.irStdev;
%% minimize least squares for 45min
a0=[2 1 -10];
lb = [0 0 -50];
ub=[5 5 0];
gDelta=9;
iDelta=9;
tic
[GlobmodelFits45min, Globstats45min]=GlobalRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc
% minimize least squares for 60min
a0=[2 1 -10];
lb = [0 0 -50];
ub=[4 4 0];
gDelta=12;
iDelta=12;
tic
[GlobmodelFits60min, Globstats60min]=GlobalRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc
% 120min
a0=[2 1 -10];
lb = [0 0 -50];
ub=[4 4 0];
gDelta=24;
iDelta=24;
tic
[GlobmodelFits120min, Globstats120min]=GlobalRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX);
toc

%% Plots
ID=1:MAX;
gDelta=6;
iDelta=6;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    if GlobmodelFits30min(i).RES >0 && ~isnan(GlobmodelFits30min(i).RES)
        GlobmodelFits30min(i).predict=[patient(i).gCGM(1:(delta+gDelta-1))' (GlobmodelFits30min(i).Fits(1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            GlobmodelFits30min(i).Fits(2)*patient(i).gCGM(delta:(end-gDelta))+ ...
            GlobmodelFits30min(i).Fits(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))']';
        
        GlobmodelFits30min(i).indivRes=abs(GlobmodelFits30min(i).predict-patient(i).gCGM);
    else
        GlobmodelFits30min(i).indivRes=0;
    end
end
%%
gDelta=9;
iDelta=9;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    if GlobmodelFits45min(i).RES >0 && ~isnan(GlobmodelFits45min(i).RES)
    GlobmodelFits45min(i).predict=[patient(i).gCGM(1:(delta+gDelta-1))' (GlobmodelFits45min(i).Fits(1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            GlobmodelFits45min(i).Fits(2)*patient(i).gCGM(delta:(end-gDelta))+ ... 
            GlobmodelFits45min(i).Fits(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))']';
        
                GlobmodelFits45min(i).indivRes=abs(GlobmodelFits45min(i).predict-patient(i).gCGM);

    else
                GlobmodelFits45min(i).indivRes=0;

    end
end

gDelta=12;
iDelta=12;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    if GlobmodelFits60min(i).RES >0 && ~isnan(GlobmodelFits60min(i).RES)
    GlobmodelFits60min(i).predict=[patient(i).gCGM(1:(delta+gDelta-1))' (GlobmodelFits60min(i).Fits(1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            GlobmodelFits60min(i).Fits(2)*patient(i).gCGM(delta:(end-gDelta))+ ... 
            GlobmodelFits60min(i).Fits(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))']';
        
                GlobmodelFits60min(i).indivRes=abs(GlobmodelFits60min(i).predict-patient(i).gCGM);

    else
                GlobmodelFits60min(i).indivRes=0;

    end
end

gDelta=24;
iDelta=24;
delta=max(gDelta,iDelta)+1;

for i=ID
    tEND=length(patient(i).gtimes);
    if GlobmodelFits120min(i).RES >0 && ~isnan(GlobmodelFits120min(i).RES)
    GlobmodelFits120min(i).predict=[patient(i).gCGM(1:(delta+gDelta-1))' (GlobmodelFits120min(i).Fits(1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            GlobmodelFits120min(i).Fits(2)*patient(i).gCGM(delta:(end-gDelta))+ ... 
            GlobmodelFits120min(i).Fits(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))']';
        
                GlobmodelFits120min(i).indivRes=abs(GlobmodelFits120min(i).predict-patient(i).gCGM);

    else
        GlobmodelFits120min(i).indivRes=0;
    end
end
%%
for i=1:20
    if ~isempty(GlobmodelFits120min(i).predict)
    figure(i)
    subplot(2,1,1)
    plot(patient(i).gtimes,patient(i).gCGM,'-o',patient(i).gtimes,GlobmodelFits30min(i).predict,patient(i).gtimes,GlobmodelFits45min(i).predict,patient(i).gtimes,GlobmodelFits60min(i).predict,patient(i).gtimes,GlobmodelFits120min(i).predict)
    legend('CGM data','30min','45min','60min','120min')
    title('CGM actual vs Predicted values for various time steps')
    ylabel('Glucose Measurements (mg/dL)')
    xlabel('minutes since start of trial')
    subplot(2,1,2)
    plot(patient(i).times,patient(i).IOB)
    ylabel('Insulin On Board')
    xlabel('minutes since start of trial')
    title(allfiles(i).name)
    else
        continue
    end
end
% %%
% for i=ID
% figure(100+i)
% subplot(2,1,1)
% hist(modelFits30min(i).predict(13:end)-patient(i).gCGM(13:end),200)
% title('Residuals 30min')
% subplot(2,1,2)
% hist(GlobmodelFits30min(i).predict(13:end)-patient(i).gCGM(13:end),200)
% title('Residuals 30min Global')
% end
%% calculate stats for residuals of global fit model

Globstats30min.indivRes=vertcat(GlobmodelFits30min(1:end).indivRes);
Globstats30min.irMean=nanmean(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.irStdev=nanstd(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.irMin=min(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.irMax=max(Globstats30min.indivRes(Globstats30min.indivRes>0));
Globstats30min.ir95=1.96*Globstats30min.irStdev;

Globstats45min.indivRes=vertcat(GlobmodelFits45min(1:end).indivRes);
Globstats45min.irMean=nanmean(Globstats45min.indivRes(Globstats45min.indivRes>0));
Globstats45min.irStdev=nanstd(Globstats45min.indivRes(Globstats45min.indivRes>0));
Globstats45min.irMin=min(Globstats45min.indivRes(Globstats45min.indivRes>0));
Globstats45min.irMax=max(Globstats45min.indivRes(Globstats45min.indivRes>0));
Globstats45min.ir95=1.96*Globstats45min.irStdev;

Globstats60min.indivRes=vertcat(GlobmodelFits60min(1:end).indivRes);
Globstats60min.irMean=nanmean(Globstats60min.indivRes(Globstats60min.indivRes>0));
Globstats60min.irStdev=nanstd(Globstats60min.indivRes(Globstats60min.indivRes>0));
Globstats60min.irMin=min(Globstats60min.indivRes(Globstats60min.indivRes>0));
Globstats60min.irMax=max(Globstats60min.indivRes(Globstats60min.indivRes>0));
Globstats60min.ir95=1.96*Globstats60min.irStdev;

Globstats120min.indivRes=vertcat(GlobmodelFits120min(1:end).indivRes);
Globstats120min.irMean=nanmean(Globstats120min.indivRes(Globstats120min.indivRes>0));
Globstats120min.irStdev=nanstd(Globstats120min.indivRes(Globstats120min.indivRes>0));
Globstats120min.irMin=min(Globstats120min.indivRes(Globstats120min.indivRes>0));
Globstats120min.irMax=max(Globstats120min.indivRes(Globstats120min.indivRes>0));
Globstats120min.ir95=1.96*Globstats120min.irStdev;

%% Optimization by window
tic 

a0=[.5 .5 -2.5];
lb = [0 0 -10];
ub=[4 4 0];
gDelta=6;
iDelta=6;
stepsz=20;
ovlp=5;

[modelFits30Win20, stats30Win20] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=40;
ovlp=5;
[modelFits30Win40, stats30Win40] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

stepsz=60;
ovlp=5;
[modelFits30Win60, stats30Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


% 45min

a0=[2 1 -10];
lb = [0 0 -50];
ub=[2 2 0];
gDelta=9;
iDelta=9;

stepsz=20;
ovlp=5;
[modelFits45Win20, stats45Win20] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=40;
ovlp=5;
[modelFits45Win40, stats45Win40] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits45Win60, stats45Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


% 60min

a0=[2 1 -10];
lb = [0 0 -50];
ub=[4 4 0];
gDelta=12;
iDelta=12;

stepsz=20;
ovlp=5;
[modelFits60Win20, stats60Win20] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=40;
ovlp=5;
[modelFits60Win40, stats60Win40] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);


stepsz=60;
ovlp=5;
[modelFits60Win60, stats60Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

%120min

a0=[2 1 -10];
lb = [0 0 -50];
ub=[4 4 0];
gDelta=24;
iDelta=24;
stepsz=20;
ovlp=5;

[modelFits120Win20, stats120Win20] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

stepsz=40;
ovlp=5;
[modelFits120Win40, stats120Win40] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);

stepsz=60;
ovlp=5;

[modelFits120Win60, stats120Win60] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp);
toc

%% sensitivity analysis
i=167;
%construct object for numeric parameter that can take specified values in
%distribution
a1 = param.Continuous('a1',stats45Win20.mean(1));
a2 = param.Continuous('a2',stats45Win20.mean(2));
a3 = param.Continuous('a3',abs(stats45Win20.mean(3)));

%make normal distrib using parameters
pdR1 = makedist('Normal','mu',stats45Win20.mean(1),'sigma',stats45Win20.stdev(1));
pdR2 = makedist('Normal','mu',stats45Win20.mean(2),'sigma',stats45Win20.stdev(2));
pdR3 = makedist('Normal','mu',stats45Win20.mean(3),'sigma',stats45Win20.stdev(3));

x1=linspace(stats45Win20.mean(1)-3*stats45Win20.stdev(1),stats45Win20.mean(1)+3*stats45Win20.stdev(1));
x2=linspace(stats45Win20.mean(2)-3*stats45Win20.stdev(2),stats45Win20.mean(2)+3*stats45Win20.stdev(2));
x3=linspace(stats45Win20.mean(3)-3*stats45Win20.stdev(3),stats45Win20.mean(3)+3*stats45Win20.stdev(3));

figure(1)
subplot(3,1,1)
plot(x1,pdf(pdR1,x1));
subplot(3,1,2)
plot(x2,pdf(pdR2,x2));
subplot(3,1,3)
plot(x3,pdf(pdR3,x3));
%all this shit to add title to the top of subplots
set(gcf,'NextPlot','add');
axes;
h = title('Parameter distributions');
set(gca,'Visible','off');
set(h,'Visible','on');

%specify pdR as prob dist for a1 param in sdo.paramspace
ps1=sdo.ParameterSpace(a1,pdR1);
ps2=sdo.ParameterSpace(a2,pdR2);
ps3=sdo.ParameterSpace(a3,pdR3);

%generate samples
Ns=20;
X1=sdo.sample(ps1,Ns);
X2=sdo.sample(ps2,Ns);
X3=sdo.sample(ps3,Ns);


%plot samples
figure(2)
subplot(2,2,1)
sdo.scatterPlot(X1,X2)
subplot(2,2,2)
sdo.scatterPlot(X2,X3)
subplot(2,2,3)
sdo.scatterPlot(X1,X3)
%title
set(gcf,'NextPlot','add');
axes;
h = title('Parameter combinations');
set(gca,'Visible','off');
set(h,'Visible','on');

% construct array of parameter combinations
clear Sens30 dist30 vC
dist30=combvec(table2array(X1)',table2array(X2)',table2array(X3)')';
%keep only valid param combos
vC=find(dist30(:,3,:)<0 & dist30(:,1,:)>0 & dist30(:,2,:)>0); %& dist30(:,1,:)+dist30(:,2,:)<=1.2 & dist30(:,1,:)+dist30(:,2,:)>=1

ID=1:MAX;
gDelta=6;
iDelta=6;
delta=max(gDelta,iDelta)+1;

Sens30(length(vC)).FitsRES=[];
Sens30(length(vC)).ogFitsRES=[];


for n=1:length(vC)
        Sens30(n).FitsRES=(patient(i).gCGM(delta+gDelta:end)-(dist30(vC(n),1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            dist30(vC(n),2)*patient(i).gCGM(delta:(end-gDelta))+ ...
            dist30(vC(n),3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta))));
        
        Sens30(n).ogFitsRES=(patient(i).gCGM(delta+gDelta:end)-((stats45Win20.mean(1))*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            stats45Win20.mean(2)*patient(i).gCGM(delta:(end-gDelta))+ ...
            stats45Win20.mean(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta))));
        
        Sens30(n).cinRES=Sens30(n).ogFitsRES-Sens30(n).FitsRES;
        Sens30(n).totcinRES=sum(abs(Sens30(n).cinRES));
%         Sens30(n).RESES=abs(Sens30(n).Fits-patient(i).gCGM);
%         Sens30(n).RES=sum(Sens30(n).RESES);
%         Sens30(n).rMean=mean(Sens30(n).RESES);
%         Sens30(n).stdev=std(Sens30(n).RESES);
%         Sens30(n).max=max(Sens30(n).RESES);
%         Sens30(n).min=min(Sens30(n).RESES(delta+gDelta:end));
end

%sort array by res mean
clear res ind porder ord percentChange
[res, ind]=sort([Sens30.totcinRES]);
porder=vC(ind); %get ordering for which param was used
ord=dist30(porder,1:3); %a(1) is most sensitive on first-pass

percentChange=[ord./repmat(stats45Win20.mean,n,1)-ones(n,1), zeros(n,1), res'];

%add weights to dots at vC based on the total residue
figure(203)
scatter3(dist30(vC,1),dist30(vC,2),dist30(vC,3),.01.*[Sens30.totcinRES],.01.*[Sens30.totcinRES],'*')
xlabel('a(1)')
ylabel('a(2)')
zlabel('a(3)')
title('Change in residue for various parameter combinations around the mean')
colorbar

%backparse to get how much change in each param affects output

%%
for i=1:20
    if length(modelFits45Win20(i).predict)==length(patient(i).gtimes(19:end))
    % plot
    figure(i)
    subplot(2,1,1)
    plot(patient(i).gtimes(19:end),patient(i).gCGM(19:end),'-o',patient(i).gtimes(19:end),modelFits45Win20(i).predict)
    legend('CGM data','45Win20 predict')
    title('CGM actual vs Predicted values for various time steps')
    ylabel('Glucose Measurements (mg/dL)')
    xlabel('minutes since start of trial')
    subplot(2,1,2)
    plot(patient(i).times,patient(i).IOB)
    ylabel('Insulin On Board')
    xlabel('minutes since start of trial')
    title(allfiles(i).name)
    end
end