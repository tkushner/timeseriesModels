function [] = writetoFileTex(fitID,pID,stats)

fileID=fopen(strcat(pID,'_clean_windowFits.txt'),'a');
fprintf(fileID,'%12s\n',fitID);
fprintf(fileID,'%8s %6s %6s %6s\n','parameter &', 'a1 &', 'a2 &', 'a3');
fprintf(fileID,'%8s %.4f %3s %.4f %3s %.4f\n','mean & ',stats.mean(1),' & ',stats.mean(2),' & ',stats.mean(3));
fprintf(fileID,'%8s %.4f %3s %.4f %3s %.4f\n','std dev & ',stats.stdev(1),' & ',stats.stdev(2),' & ',stats.stdev(3));
fprintf(fileID, '%30s\n','Residual Mean & Stdev & 95\% Confidence & Min & Max');
fprintf(fileID,'%.4f %3s %.4f %6s %.4f %3s %.4f %3s %.4f\n',stats.RESmean, ' & ',stats.RESstdev, '& $\pm$ ', stats.RES95, ' & ', stats.AbsMIN, ' & ', stats.AbsMAX);
fclose(fileID);
end
