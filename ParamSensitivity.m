function [out] = ParamSensitivity(in)

%% sensitivity analysis
i=1;
%construct object for numeric parameter that can take specified values in
%distribution
a1 = param.Continuous('a1',GlobmodelFits30min(i).Fits(1));
a2 = param.Continuous('a2',GlobmodelFits30min(i).Fits(2));
a3 = param.Continuous('a3',abs(GlobmodelFits30min(i).Fits(3)));

%create pdR to configure param space
x1 = [0.95 0.99 1.01 1.05]*a1.Value;
x2 = [0.95 0.99 1.01 1.05]*a2.Value;   
x3 = [0.95 0.99 1.01 1.05]*a3.Value;  
F = [0 0.5 0.5 1];

%make normal distrib using parameters
pdR1 = makedist('Normal','mu',Globstats30min.mean(1),'sigma',Globstats30min.stdev(1));
pdR2 = makedist('Normal','mu',Globstats30min.mean(2),'sigma',Globstats30min.stdev(2));
pdR3 = makedist('Normal','mu',Globstats30min.mean(3),'sigma',Globstats30min.stdev(3));


x1=linspace(Globstats30min.mean(1)-3*Globstats30min.stdev(1),Globstats30min.mean(1)+3*Globstats30min.stdev(1));
x2=linspace(Globstats30min.mean(2)-3*Globstats30min.stdev(2),Globstats30min.mean(2)+3*Globstats30min.stdev(2));
x3=linspace(Globstats30min.mean(3)-3*Globstats30min.stdev(3),Globstats30min.mean(3)+3*Globstats30min.stdev(3));

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

%% construct array of parameter combinations
dist30=combvec(table2array(X1)',table2array(X2)',table2array(X3)')';
%keep only valid param combos
vC=find(dist30(:,3,:)<0 & dist30(:,1,:)>0 & dist30(:,2,:)>0);

gDelta=6;
iDelta=6;
delta=max(gDelta,iDelta)+1;

Sens30(length(vC)).Fits=[];
Sens30(length(vC)).RESES=[];
Sens30(length(vC)).RES=[];

for n=1:length(vC)
        Sens30(n).Fits=[patient(i).gCGM(1:(delta+gDelta-1))' (dist30(vC(n),1)*patient(i).gCGM((delta-gDelta):(end-2*gDelta))+ ...
            dist30(vC(n),2)*patient(i).gCGM(delta:(end-gDelta))+ ...
            dist30(vC(n),3)*patient(i).gIOB((delta-iDelta):(end-2*iDelta)))']';
        
        Sens30(n).RESES=abs(Sens30(n).Fits-patient(i).gCGM);
        Sens30(n).RES=sum(Sens30(n).RESES);
        Sens30(n).rMean=mean(Sens30(n).RESES);
        Sens30(n).stdev=std(Sens30(n).RESES);
        Sens30(n).max=max(Sens30(n).RESES);
        Sens30(n).min=min(Sens30(n).RESES(delta+gDelta:end));
end
%%
%add weights to dots at vC based on the total residue
figure(203)
scatter3(dist30(vC,1),dist30(vC,2),dist30(vC,3),.01.*[Sens30.RES],.01.*[Sens30.RES],'*')
xlabel('a(1)')
ylabel('a(2)')
zlabel('a(3)')



end
