function [modelFits, stats] = RegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX)
delta=max(gDelta,iDelta)+1;

%initialize struct for modelFit
modelFits(MAX).Fits=[];
modelFits(MAX).RES=[];
modelFits(MAX).mean=[];
modelFits(MAX).stdev=[];
modelFits(MAX).RESmean=[];
modelFits(MAX).RESstdev=[];
modelFits(MAX).RESmin=[];

for i = 1:MAX
        tEND=length(patient(i).gCGM)-delta;
        modelFits(i).Fits=zeros(tEND,3);
        modelFits(i).RES=zeros(tEND,1);
end

parfor i = 1:MAX
    
    tEND=length(patient(i).gCGM)-delta;
    %check for trials with no errors (drop out/ extra points / etc)
    if min(diff(patient(i).gtimes))==5 && max(diff(patient(i).gtimes)==5)
        n=1;
        for t=delta:tEND-1
            fun=@(a)(patient(i).gCGM(t+gDelta)-a(1)*patient(i).gCGM(t-gDelta)-a(2)*patient(i).gCGM(t)-a(3)*patient(i).gIOB(t-iDelta))^2;
            [a,Res]=fmincon(fun,a0,[],[],[],[],lb,ub);
            modelFits(i).Fits(n,:)=a;
            modelFits(i).RES(n)=sqrt(Res);
            n=n+1;
        end
        %Find first zero entry (to overcome fact that you cannot initialize empty matrix)
        tENDfits=min(find(modelFits(i).Fits==[0,0,0]))-1;
        %look at means & std dev for each parameter
        modelFits(i).mean=mean(modelFits(i).Fits(1:tENDfits,1:3));
        modelFits(i).stdev=std(modelFits(i).Fits(1:tENDfits,1:3));
        modelFits(i).RESmean=mean(modelFits(i).RES(1:tENDfits));
        modelFits(i).RESstdev=std(modelFits(i).RES(1:tENDfits));
        modelFits(i).RESmin=min(modelFits(i).RES(1:tENDfits));

    elseif sum(diff(patient(i).gtimes)~=5)==1
        n=1;
        SWITCH=find(diff(patient(i).gtimes)~=5);
        tEND1=SWITCH-delta;
        tEND2=length(patient(i).gCGM)-delta;
        if (tEND1-1)>delta
            for t=delta:(tEND1-1)
            fun=@(a)(patient(i).gCGM(t+gDelta)-a(1)*patient(i).gCGM(t-gDelta)-a(2)*patient(i).gCGM(t)-a(3)*patient(i).gIOB(t-iDelta))^2;
                [a,Res]=fmincon(fun,a0,[],[],[],[],lb,ub);
                modelFits(i).Fits(n,:)=a;
                modelFits(i).RES(n)=sqrt(Res);
                n=n+1;
            end
            for t=(SWITCH+delta+1):(tEND2-1)
            fun=@(a)(patient(i).gCGM(t+gDelta)-a(1)*patient(i).gCGM(t-gDelta)-a(2)*patient(i).gCGM(t)-a(3)*patient(i).gIOB(t-iDelta))^2;
                [a,Res]=fmincon(fun,a0,[],[],[],[],lb,ub);
                modelFits(i).Fits(n,:)=a;
                modelFits(i).RES(n)=sqrt(Res);
                n=n+1;
            end
        else
            for t = (SWITCH+delta+1):(tEND2-1)
            fun=@(a)(patient(i).gCGM(t+gDelta)-a(1)*patient(i).gCGM(t-gDelta)-a(2)*patient(i).gCGM(t)-a(3)*patient(i).gIOB(t-iDelta))^2;
                [a,Res]=fmincon(fun,a0,[],[],[],[],lb,ub);
                modelFits(i).Fits(n,:)=a;
                modelFits(i).RES(n)=sqrt(Res);
                n=n+1;
            end
        end
        %Find first zero entry (to overcome fact that you cannot initialize empty matrix)
        tENDfits=min(find(modelFits(i).Fits==[0,0,0]))-1;
        %look at means & std dev for each parameter
        modelFits(i).mean=mean(modelFits(i).Fits(1:tENDfits,1:3));
        modelFits(i).stdev=std(modelFits(i).Fits(1:tENDfits,1:3));
        modelFits(i).RESmean=mean(modelFits(i).RES(1:tENDfits));
        modelFits(i).RESstdev=std(modelFits(i).RES(1:tENDfits));
        modelFits(i).RESmin=min(modelFits(i).RES(1:tENDfits));
    else continue
    end
    %check if empty then set as NaN (else issues with padcat)
    if isempty(modelFits(i).RESmin)
        modelFits(i).RESmin=NaN(1)
    end
end
[MEAN]=padcat(modelFits(1:end).mean);
[RESMEAN]=padcat(modelFits(1:end).RESmean);
stats.mean=nanmean(MEAN);
stats.stdev=nanstd(MEAN);
stats.RESmean=nanmean(RESMEAN);
stats.RESstdev=nanstd(RESMEAN);
stats.RESmax=max(max(padcat(modelFits(1:end).RES)));
stats.RESmin=min(min(padcat(modelFits(1:end).RESmin)));
stats.RES95=1.96*stats.RESstdev;
end