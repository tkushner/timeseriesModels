%clustering
%Taisa Kushner
%1 Aug 2017

function [] = clusterFits(dataset, datadef, numclust, plotsOn)

%specifiy dataset and plot label
%dataset=modelFits30Win36;
%datadef='Delta=30min, fit window=180min';

allfitsnan=vertcat(vertcat(dataset.Fits));
windS=horzcat(dataset.windowS);
windE=horzcat(dataset.windowE);
allfitsTimeS=datetime(windS(~cellfun('isempty',windS)),'InputFormat','HH:mm:ss MM/dd/yyyy ');
allfitsTimeE=datetime(windE(~cellfun('isempty',windE)),'InputFormat','HH:mm:ss MM/dd/yyyy ');


allfits=~isnan(allfitsnan);
sSize=allfitsTimeS.Hour+allfitsTimeS.Minute/60;
eSize=allfitsTimeE.Hour+allfitsTimeE.Minute/60;
morn=find(sSize<=13);
morn2=find(eSize<=13);
sSize(morn)=sSize(morn)+24;
eSize(morn2)=eSize(morn2)+24;

A=allfitsnan(allfits(:,1),1);
B=allfitsnan(allfits(:,2),2);
C=allfitsnan(allfits(:,3),3);

if numclust==3
    Xa1=[sSize',A];
    
    opts = statset('Display','final');
    [idxa1,Csa1] = kmeans(Xa1,numclust,'Distance','cityblock',...
        'Replicates',8,'Options',opts);
    
    figure(102)
    subplot(2,1,1)
    plot(Xa1(idxa1==1,1),Xa1(idxa1==1,2),'r.','MarkerSize',12)
    hold on
    plot(Xa1(idxa1==2,1),Xa1(idxa1==2,2),'b.','MarkerSize',12)
    plot(Xa1(idxa1==3,1),Xa1(idxa1==3,2),'g.','MarkerSize',12)
    plot(Csa1(:,1),Csa1(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3)
    legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW')
    %legend('Cluster 1','Cluster 2','Centroids','Location','NW')
    title(strcat({'Cluster Assignments and Centroids for '},datadef))
    xlabel('time, hrs')
    ylabel('a(1)')
    hold off
    %
    Xa2=[sSize',B];
    
    opts = statset('Display','final');
    [idxa2,Csa2] = kmeans(Xa2,numclust,'Distance','cityblock',...
        'Replicates',8,'Options',opts);
    
    figure(102)
    subplot(2,1,2)
    plot(Xa2(idxa2==1,1),Xa2(idxa2==1,2),'r.','MarkerSize',12)
    hold on
    plot(Xa2(idxa2==2,1),Xa2(idxa2==2,2),'b.','MarkerSize',12)
    plot(Xa2(idxa2==3,1),Xa2(idxa2==3,2),'g.','MarkerSize',12)
    plot(Csa2(:,1),Csa2(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3)
    legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW')
    %legend('Cluster 1','Cluster 2','Centroids','Location','NW')
    title(strcat({'Cluster Assignments and Centroids for '},datadef))
    xlabel('time, hrs')
    ylabel('a(2)')
    hold off
    % find distributions
    a1clust1=fitdist(Xa1(idxa1==1,2),'Normal');
    a1clust2=fitdist(Xa1(idxa1==2,2),'Normal');
    a1clust3=fitdist(Xa1(idxa1==3,2),'Normal');
    
    a2clust1=fitdist(Xa2(idxa2==1,2),'Normal');
    a2clust2=fitdist(Xa2(idxa2==2,2),'Normal');
    a2clust3=fitdist(Xa2(idxa2==3,2),'Normal');
    
    %determine if two datasets come from the same distribution - based on the clustering between a1 and time -- if true, null
    %hypothesis is rejected - ksp will print what the actual alpha is
    [a13ksh, a13ksp]=kstest2(Xa1(idxa1==1,2),Xa1(idxa1==3,2),'Alpha',.35);
    [a12ksh, a12ksp]=kstest2(Xa1(idxa1==1,2),Xa1(idxa1==2,2),'Alpha',.35);
    [a23ksh, a23ksp]=kstest2(Xa1(idxa1==2,2),Xa1(idxa1==3,2),'Alpha',.35);
    
    [a13ksh2, a13ksp2]=kstest2(Xa2(idxa2==1,2),Xa2(idxa2==3,2),'Alpha',.35);
    [a12ksh2, a12ksp2]=kstest2(Xa2(idxa2==1,2),Xa2(idxa2==2,2),'Alpha',.35);
    [a23ksh2, a23ksp2]=kstest2(Xa2(idxa2==2,2),Xa2(idxa2==3,2),'Alpha',.35);
    
elseif numclust==2
    Xa1=[sSize',A];
    
    opts = statset('Display','final');
    [idxa1,Csa1] = kmeans(Xa1,numclust,'Distance','cityblock',...
        'Replicates',8,'Options',opts);
    
    figure(102)
    subplot(2,1,1)
    plot(Xa1(idxa1==1,1),Xa1(idxa1==1,2),'r.','MarkerSize',12)
    hold on
    plot(Xa1(idxa1==2,1),Xa1(idxa1==2,2),'b.','MarkerSize',12)
    plot(Csa1(:,1),Csa1(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3)
    legend('Cluster 1','Cluster 2','Centroids','Location','NW')
    title(strcat({'Cluster Assignments and Centroids for '},datadef))
    xlabel('time, hrs')
    ylabel('a(1)')
    hold off
    %
    Xa2=[sSize',B];
    
    opts = statset('Display','final');
    [idxa2,Csa2] = kmeans(Xa2,numclust,'Distance','cityblock',...
        'Replicates',8,'Options',opts);
    
    figure(102)
    subplot(2,1,2)
    plot(Xa2(idxa2==1,1),Xa2(idxa2==1,2),'r.','MarkerSize',12)
    hold on
    plot(Xa2(idxa2==2,1),Xa2(idxa2==2,2),'b.','MarkerSize',12)
    plot(Xa2(idxa2==3,1),Xa2(idxa2==3,2),'g.','MarkerSize',12)
    plot(Csa2(:,1),Csa2(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3)
    legend('Cluster 1','Cluster 2','Centroids','Location','NW')
    title(strcat({'Cluster Assignments and Centroids for '},datadef))
    xlabel('time, hrs')
    ylabel('a(2)')
    hold off
    % find distributions
    a1clust1=fitdist(Xa1(idxa1==1,2),'Normal');
    a1clust2=fitdist(Xa1(idxa1==2,2),'Normal');
    
    a2clust1=fitdist(Xa2(idxa2==1,2),'Normal');
    a2clust2=fitdist(Xa2(idxa2==2,2),'Normal');
    
    %determine if two datasets come from the same distribution - based on the clustering between a1 and time -- if true, null
    %hypothesis is rejected - ksp will print what the actual alpha is
    [a12ksh, a12ksp]=kstest2(Xa1(idxa1==1,2),Xa1(idxa1==2,2),'Alpha',.35);
    
    [a12ksh2, a12ksp2]=kstest2(Xa2(idxa2==1,2),Xa2(idxa2==2,2),'Alpha',.35);
    
end


%plot clusters in 3d
INDEX3D=40*ones(length(A),1);
INDEX3D(find(Xa1(idxa1==1,1)))=10;
INDEX3D(find(Xa1(idxa1==2,1)))=25;

if plotsOn==1
    figure(103)
    scatter3(A,B,C,sSize,INDEX3D)
    xlabel('a(1)')
    ylabel('a(2)')
    zlabel('a(3)')
    title(strcat({'Parameter fits for '}, datadef, {' clustered by a1 parameter in time'}))
    
    
    %method 2 for clustering
    clear AA
    AA=clusterdata([A,B,C],'criterion','distance','linkage','ward','maxclust',3); %ward is inner squared distance
    figure(104)
    scatter3(A,B,C,INDEX3D,AA,'filled')
    xlabel('a(1)')
    ylabel('a(2)')
    zlabel('a(3)')
    title(strcat('Clustering of all parameter fits for .', datadef))
    
    figure(105)
    subplot(3,3,1)
    scatter3(A,B,sSize)
    xlabel('a(1) parameter')
    ylabel('a(2) parameter')
    zlabel('start time for window')
    subplot(3,3,2)
    plot(sSize,A,'o')
    xlabel('start time for window')
    ylabel('a(1) parameter')
    subplot(3,3,3)
    plot(sSize,B,'o')
    xlabel('start time for window')
    ylabel('a(2) parameter')
    subplot(3,3,4)
    scatter3(A,B,eSize)
    xlabel('a(1) parameter')
    ylabel('a(2) parameter')
    zlabel('end time for window')
    subplot(3,3,5)
    plot(eSize,A,'o')
    xlabel('end time for window')
    ylabel('a(1) parameter')
    subplot(3,3,6)
    plot(eSize,B,'o')
    xlabel('end time for window')
    ylabel('a(2) parameter')
    subplot(3,3,7)
    plot(sSize,C,'o')
    xlabel('start time for window')
    ylabel('a(3) parameter')
    subplot(3,3,8)
    plot(eSize,C,'o')
    xlabel('end time for window')
    ylabel('a(3) parameter')
    %title
    set(gcf,'NextPlot','add');
    axes;
    h = title(datadef);
    set(gca,'Visible','off');
    set(h,'Visible','on');
end

end

