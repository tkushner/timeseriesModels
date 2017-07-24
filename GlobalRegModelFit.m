%optimize across the entire function rather than at each time interval
function [modelFits, stats] = GlobalRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX)

delta=max(gDelta,iDelta)+1;

%initialize struct for modelFit
modelFits(MAX).Fits=[];
modelFits(MAX).RES=[];

for i = 1:MAX
    tBEG=delta+gDelta;
    modelFits(i).Fits=NaN(1,3);
    modelFits(i).RES=NaN(1);
    %check for trials with no errors (drop out/ extra points / etc)
    if max(diff(patient(i).gtimes))==5 & min(diff(patient(i).gtimes)==5)
        fun=@(a)(sum((patient(i).gCGM(tBEG:end)-(a(1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            a(2)*patient(i).gCGM(delta:(end-gDelta))+a(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))).^2));
         options=optimoptions('Display','off');
        [a,Res]=fmincon(fun,a0,[],[],[],[],lb,ub);
        modelFits(i).Fits=a;
        modelFits(i).RES=sqrt(Res)\length(patient(i).gCGM(tBEG:end));          
    else
        continue
    end
    %look at means & std dev for each parameter
end
[MEAN]=padcat(modelFits(1:end).Fits);
[RESES]=padcat(modelFits(1:end).RES);
stats.mean=nanmean(MEAN);
stats.stdev=nanstd(MEAN);
stats.RESmean=nanmean(RESES);
stats.RESstdev=nanstd(RESES);
stats.RESmax=max(RESES);
stats.RESmin=min(RESES);
stats.RES95=1.96*stats.RESstdev;
end
