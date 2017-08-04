%calculate residuals and other data for cluster fits
%Taisa Kushner
%4 August 2017

function [clustData, clustStats]=calcClustData(Early, Mid, Late, gDelta, iDelta, patient,MAX,M,L,numClust);
delta=max(gDelta,iDelta)+1;

%separate data into early/mid/late
clear clustData clustStats
clustData(MAX).resEarly=[];
clustData(MAX).resMid=[];
clustData(MAX).resLate=[];

if numClust==3
    for i=1:MAX
        clear cgmEarly cgmMid cgmLate iobEarly iobMid iobLate hours tmstmp
        tmstmp=datetime(patient(i).gDatetime, 'InputFormat','HH:mm:ss MM/dd/yyyy ');
        hours=tmstmp.Hour;
        hours(find(tmstmp.Hour<13))=hours(find(tmstmp.Hour<13))+24;
        cgmEarly=patient(i).gCGM(find(hours<M));
        cgmMid=patient(i).gCGM(find(hours>=M & hours<L));
        cgmLate=patient(i).gCGM(find(hours>=L));
        iobEarly=patient(i).gIOB(find(hours<M));
        iobMid=patient(i).gIOB(find(hours>=M & hours<L));
        iobLate=patient(i).gIOB(find(hours>=L));
        
        clustData(i).resEarly=cgmEarly(delta+gDelta:end)-((Early(1))*cgmEarly((delta-gDelta):(end-2*gDelta))+ ...
            Early(2)*cgmEarly(delta:(end-gDelta))+ ...
            Early(3)*iobEarly((delta-iDelta):(end-2*iDelta)));
        
        clustData(i).resMid=cgmMid(delta+gDelta:end)-((Mid(1))*cgmMid((delta-gDelta):(end-2*gDelta))+ ...
            Mid(2)*cgmMid(delta:(end-gDelta))+ ...
            Mid(3)*iobMid((delta-iDelta):(end-2*iDelta)));
        
        
        clustData(i).resLate=cgmLate(delta+gDelta:end)-((Late(1))*cgmLate((delta-gDelta):(end-2*gDelta))+ ...
            Late(2)*cgmLate(delta:(end-gDelta))+ ...
            Late(3)*iobLate((delta-iDelta):(end-2*iDelta)));
        
    end
    clustStats.meanEarly=nanmean(vertcat(clustData.resEarly));
    clustStats.meanMid=nanmean(vertcat(clustData.resMid));
    clustStats.meanLate=nanmean(vertcat(clustData.resLate));
    
    clustStats.stdEarly=nanstd(vertcat(clustData.resEarly));
    clustStats.stdMid=nanstd(vertcat(clustData.resMid));
    clustStats.stdLate=nanstd(vertcat(clustData.resLate));
    
    clustStats.early95=1.96*clustStats.stdEarly;
    clustStats.mid95=1.96*clustStats.stdMid;
    clustStats.late95=1.96*clustStats.stdLate;
    
    clustStats.ABSmeanEarly=nanmean(abs(vertcat(clustData.resEarly)));
    clustStats.ABSmeanMid=nanmean(abs(vertcat(clustData.resMid)));
    clustStats.ABSmeanLate=nanmean(abs(vertcat(clustData.resLate)));
    
    clustStats.ABSstdEarly=nanstd(abs(vertcat(clustData.resEarly)));
    clustStats.ABSstdMid=nanstd(abs(vertcat(clustData.resMid)));
    clustStats.ABSstdLate=nanstd(abs(vertcat(clustData.resLate)));
    
    clustStats.ABSearly95=1.96*clustStats.ABSstdEarly;
    clustStats.ABSmid95=1.96*clustStats.ABSstdMid;
    clustStats.ABSlate95=1.96*clustStats.ABSstdLate;
    
elseif numClust==2
    
    for i=1:MAX
        clear cgmEarly cgmMid cgmLate iobEarly iobMid iobLate hours tmstmp
        tmstmp=datetime(patient(i).gDatetime, 'InputFormat','HH:mm:ss MM/dd/yyyy ');
        hours=tmstmp.Hour;
        hours(find(tmstmp.Hour<13))=hours(find(tmstmp.Hour<13))+24;
        cgmEarly=patient(i).gCGM(find(hours<M));
        cgmLate=patient(i).gCGM(find(hours>=L));
        iobEarly=patient(i).gIOB(find(hours<M));
        iobLate=patient(i).gIOB(find(hours>=L));
        
        clustData(i).resEarly=cgmEarly(delta+gDelta:end)-((Early(1))*cgmEarly((delta-gDelta):(end-2*gDelta))+ ...
            Early(2)*cgmEarly(delta:(end-gDelta))+ ...
            Early(3)*iobEarly((delta-iDelta):(end-2*iDelta)));
        
        clustData(i).resLate=cgmLate(delta+gDelta:end)-((Late(1))*cgmLate((delta-gDelta):(end-2*gDelta))+ ...
            Late(2)*cgmLate(delta:(end-gDelta))+ ...
            Late(3)*iobLate((delta-iDelta):(end-2*iDelta)));
        
    end
    clustStats.meanEarly=nanmean(vertcat(clustData.resEarly));
    clustStats.meanLate=nanmean(vertcat(clustData.resLate));
    
    clustStats.stdEarly=nanstd(vertcat(clustData.resEarly));
    clustStats.stdLate=nanstd(vertcat(clustData.resLate));
    
    clustStats.early95=1.96*clustStats.stdEarly;
    clustStats.late95=1.96*clustStats.stdLate;
    
    clustStats.ABSmeanEarly=nanmean(abs(vertcat(clustData.resEarly)));
    clustStats.ABSmeanLate=nanmean(abs(vertcat(clustData.resLate)));
    
    clustStats.ABSstdEarly=nanstd(abs(vertcat(clustData.resEarly)));
    clustStats.ABSstdLate=nanstd(abs(vertcat(clustData.resLate)));
    
    clustStats.ABSearly95=1.96*clustStats.ABSstdEarly;
    clustStats.ABSlate95=1.96*clustStats.ABSstdLate;
end
end