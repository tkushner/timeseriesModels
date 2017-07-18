%optimize across windows within each trial
%function [modelFits, stats] = WindowRegModelFit(a0,lb,ub,gDelta,iDelta,patient,MAX)

a0=[.5 .5 -2.5];
lb = [0 0 -10];
ub=[1 1 0];
gDelta=6;
iDelta=6;
MAX=10;

delta=max(gDelta,iDelta)+1;

%initialize struct for modelFit
modelFits(MAX).Fits=[];
modelFits(MAX).RES=[];
modelFits(MAX).numwin=[];

%define windows for parsing & overlap (default: 300min (60 steps), 25min
%overlap (5)
stepsz= 60;
ovlp = 5;
step=stepsz-ovlp;

for i = 1:MAX
    clear CGM;
    clear IOB;
    
    %determine how many windows
    Nw=floor(length(patient(i).gtimes)/step);
    modelFits(i).numwin=Nw;
    modelFits(i).Fits=NaN(Nw,3);
    modelFits(i).RES=NaN(Nw,1);
    
    %check no drop out or double cgm readings
    if min(diff(patient(i).gtimes))==5 && max(diff(patient(i).gtimes)==5)
        %check at least one window
        if Nw>0
            %construct matrix of windows
            CGM=zeros(stepsz,Nw);
            IOB=zeros(stepsz,Nw);
            
            CGM(:,1)=patient(i).gCGM(1:stepsz);
            IOB(:,1)=patient(i).gIOB(1:stepsz);
            for j=1:Nw-1
                CGM(:,j+1)=patient(i).gCGM((step*j):(step-1+stepsz*j));
                IOB(:,j+1)=patient(i).gIOB((step*j):(step-1+stepsz*j));
            end
            for j=1:Nw
                fun=@(a)(sum((CGM(delta+gDelta:end,j)-(a(1)*CGM((delta-gDelta):(end-2*gDelta),j)+ ...
                    a(2)*CGM(delta:(end-gDelta),j)+a(3)*IOB((delta-iDelta):(end-2*iDelta),j))).^2));
                [a,Res]=fmincon(fun,a0,[-1,-1,0],[-1],[],[],lb,ub);
                modelFits(i).Fits(j,:)=a;
                modelFits(i).RES(j,:)=sqrt(Res)/stepsz;
            end
            
            
        
        else
            continue
        end
        
    else
        continue
    end
end
%
%     %check for trials with no errors (drop out/ extra points / etc)
%     if min(diff(patient(i).gtimes))==5 && max(diff(patient(i).gtimes)==5)
%         fun=@(a)(sum((patient(i).gCGM((delta+gDelta):end)-(a(1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
%             a(2)*patient(i).gCGM(delta:(end-gDelta))+a(3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))).^2));
%         [a,Res]=fmincon(fun,a0,[],[],[],[],lb,ub);
%         modelFits(i).Fits=a;
%         modelFits(i).RES=sqrt(Res);
%     else
%         continue
%     end
%     %look at means & std dev for each parameter
% end
% [MEAN]=padcat(modelFits(1:end).Fits);
% [RESES]=padcat(modelFits(1:end).RES);
% stats.mean=nanmean(MEAN);
% stats.stdev=nanstd(MEAN);
% stats.RESmean=nanmean(RESES);
% stats.RESstdev=nanstd(RESES);
% stats.RESmax=max(RESES);
% stats.RESmin=min(RESES);
% stats.RES95=1.96*stats.RESstdev;
%end
