%optimize across windows within each trial using nls instead of fmincon
%currently it only works for gDelta=iDelta
%Taisa Kushner 18 July 2017

function [modelFits, stats] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX,stepsz,ovlp)

% a0=[.5 .5 -2.5];
% lb = [0 0 -10];
% ub=[2 2 -1];
% gDelta=6;
% iDelta=6;
% MAX=30;

delta=max(gDelta,iDelta)+1;

%initialize struct for modelFit
modelFits(MAX).Fits=[];
modelFits(MAX).RES=[];
modelFits(MAX).numwin=[];
modelFits(MAX).ResMin=[];
modelFits(MAX).ResMax=[];
modelFits(MAX).ResMean=[];
modelFits(MAX).ResStd=[];
modelFits(MAX).MEAN=[];
modelFits(MAX).STD=[];
modelFits(MAX).predict=[];
modelFits(MAX).windowS=cell(1,1);
modelFits(MAX).windowE=cell(1,1);
modelFits(MAX).allRES=[];

%define windows for parsing & overlap (default: 300min (60 steps), 25min
%overlap (5)
% stepsz= 20;
% ovlp = 5;
step=stepsz-ovlp;

for i = 1:MAX
    clear CGM;
    clear IOB;
    
    %determine how many windows
    Nw=floor((length(patient(i).gtimes)-ovlp)/step);
    modelFits(i).numwin=Nw;
    
    %check no drop out or double cgm readings
    if (max(diff(patient(i).gtimes))==5 & min(diff(patient(i).gtimes)==5) & length(patient(i).gtimes)>stepsz)
        %check at least one window
        if Nw>0
            %construct matrix of windows
            CGM=nan(stepsz,Nw);
            IOB=nan(stepsz,Nw);
            w1=cell(Nw,1);
            w2=cell(Nw,1);
            
            CGM(:,1)=patient(i).gCGM(1:stepsz);
            IOB(:,1)=patient(i).gIOB(1:stepsz);
            w1(1)=patient(i).gDatetime(1);
            w2(1)=patient(i).gDatetime(stepsz);
            
            if Nw>1
                for j=1:Nw-1
                    CGM(:,j+1)=patient(i).gCGM((step*j):(step*j+stepsz-1));
                    IOB(:,j+1)=patient(i).gIOB((step*j):(step*j+stepsz-1));
                    w1(j+1)=patient(i).gDatetime(step*j);
                    w2(j+1)=patient(i).gDatetime(step*j+stepsz-1);
                end
            else
                continue
            end
            
            for j=1:Nw
                %construct matrix
                B=[CGM(delta+gDelta:end,j);1];
                A=[CGM(1:end-2*gDelta,j), CGM(gDelta+1:end-gDelta,j), IOB(1:end-2*iDelta,j);[1,1,0]];

                options = optimoptions('lsqlin','Algorithm','interior-point','display','off');
                [x, ~, residual, ~, ~, ~]=lsqlin(A,B,[],[],[],[],lb,ub,a0,options);
                modelFits(i).Fits(j,:)=x;
                modelFits(i).RES(j,:)=residual;
                modelFits(i).ResMin=min(modelFits(i).RES);
                modelFits(i).ResMax=max(modelFits(i).RES);
                modelFits(i).ResMean=nanmean(modelFits(i).RES);
                modelFits(i).ResStd=nanstd(modelFits(i).RES);
                modelFits(i).windowS(j)=w1(j);
                modelFits(i).windowE(j)=w2(j);
                
                if Nw>1
                    modelFits(i).MEAN=nanmean(modelFits(i).Fits);
                    modelFits(i).STD=nanstd(modelFits(i).Fits);
                else
                    modelFits(i).MEAN=nan(1,3);
                    modelFits(i).STD=NaN;
                end
            end
        else
            continue
        end
    else
        continue
    end
    
end
for i=1:length(modelFits)
    dim=size(modelFits(i).Fits);
    %check if empty
    if isempty(modelFits(i).Fits)
        modelFits(i).Fits=nan(1,3);
        modelFits(i).RES=nan(1,stepsz-2*delta+3);
    else
        for n=1:dim(1)
            A=[patient(i).gCGM(1:end-2*gDelta), patient(i).gCGM(gDelta+1:end-gDelta), patient(i).gIOB(1:end-2*iDelta)];
            x=[modelFits(i).Fits(n,1); modelFits(i).Fits(n,2); modelFits(i).Fits(n,3)];
            modelFits(i).predict(n,:)=A*x;
        end
    end
    modelFits(i).allRES=reshape(modelFits(i).RES,[1,numel(modelFits(i).RES)]);
    if isnan(modelFits(i).allRES)
        modelFits(i).allRES=[];
    end
    
end
%calculate stats
[MEAN]=vertcat(modelFits(1:end).Fits);
[RESESall]=padcat(modelFits(1:end).allRES);
RESES=reshape(RESESall,[1,numel(RESESall)]);
stats.mean=nanmean(MEAN);
stats.stdev=nanstd(MEAN);
stats.RESmean=nanmean(RESES);
stats.RESstdev=nanstd(RESES);
stats.RESmax=max(RESES);
stats.RESmin=min(RESES);
%stats.RES95=1.96*stats.RESstdev;
stats.RES95=3*stats.RESstdev;
stats.AbsMIN=min(abs(RESES));
stats.AbsMEAN=nanmean(abs(RESES));
stats.AbsMAX=max(abs(stats.RESmax),abs(stats.RESmin));
stats.AbsSTD=nanstd(abs(RESES));
stats.Abs95=1.96*stats.AbsSTD;
end


