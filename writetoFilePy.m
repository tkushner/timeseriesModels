function [] = writetoFilePy(fitID,pID,stats,delta,wind)

fileID=fopen(strcat(pID,'_pyInputs.txt'),'a');
fprintf(fileID,'%12s\n', strcat('#',fitID));
fprintf(fileID,'%8s %.4f %3s %.4f %3s %.4f %3s %d %3s %d %3s %.4f %3s %.4f %3s %.4f %3s\n', ... 
    'self.add_glucose_equation(', stats.mean(1), ',', stats.mean(2), ',', stats.mean(3), ',', ... 
    delta, ',', delta, ',', -1*stats.RES95, ',', stats.RES95, ',', stats.RESmean,')');
fclose(fileID);

%%this is for CSV
% fileID=fopen(strcat(pID,'_',fitID,'_Averaged_pyInputs.csv'),'a');
% fprintf(fileID,'%.4f %3s %.4f %3s %.4f %3s %d %3s %d %3s %.4f %3s %.4f %3s %.4f\n', ... 
%     stats.mean(1), ',', stats.mean(2), ',', stats.mean(3), ',', ... 
%     delta, ',', delta, ',', -1*stats.RES95, ',', stats.RES95, ',', stats.RESmean);
% fclose(fileID);
end
