function [] = writetoFilePyAll(fitID,pID,modelfits,stats,delta)

allfits=vertcat(modelfits.Fits);
valid=find(~isnan(allfits(:,1)));

fileID=fopen(strcat(pID,'_',fitID,'_pyInputs.csv'),'a');
for i=1:length(valid)
    fprintf(fileID,'%.4f %3s %.4f %3s %.4f %3s %d %3s %d %3s %.4f %3s %.4f %3s %.4f\n', ...
        allfits(valid(i),1), ',', allfits(valid(i),2), ',', allfits(valid(i),3), ',', ...
        delta, ',', delta, ',', -1*stats.RES95, ',', stats.RES95, ',', stats.RESmean);
end
fclose(fileID);
end
