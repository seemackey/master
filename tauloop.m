clear all;

%% Import data from spreadsheet
% using this spreadsheet:
%    Workbook: C:\Users\seema\Box\ChaseMackeyFolder\RamLab\Physiology\IC\Duration\DurationNeuro.xlsx
%    Worksheet: all


% Setup the Import Options 
opts = spreadsheetImportOptions("NumVariables", 16);

% Specify sheet and range
opts.Sheet = "all";
opts.DataRange = "A2:P1401";

% Specify column names and types
opts.VariableNames = ["Monkey", "Tank", "Block", "Loc", "comment", "cf", "cdflvl", "latency", "Threshold", "Dur", "slope", "Nthreshold", "N2threshold", "Nslope", "N2slope", "Nlatency"];
opts.VariableTypes = ["double", "string", "categorical", "categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];


% Specify variable properties
opts = setvaropts(opts, ["Tank", "comment"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Tank", "Block", "Loc", "comment"], "EmptyFieldRule", "auto");

% Import the data, need to change the path for the Mac
DurationNeuroS1 = readtable("C:\Users\seema\Box\ChaseMackeyFolder\RamLab\Physiology\IC\Duration\DurationNeuro.xlsx", opts, "UseExcel", false);

%% behavioral data

tau_beh=[26.4900
   17.1700
   30.4900
   19.1800
   18.3500
   18.0600
   17.5600
   16.6800
   19.2600
   19.3000
   17.1800
   25.9100
   17.5700
   22.3700
   22.6400
   19.8000
    9.4810
   11.8300
   26.3500
   29.4000
   18.3700
   30.3000
   26.6200
   31.3200
   32.4100
   20.2800
   14.2200
   18.4700
   16.4900
   16.3800
   19.3500
   14.7800
   33.5800
   28.8200
   17.0400
    9.6070];

% noise

tauN_beh=[18.8300
   15.4800
   14.8100
   16.9100
   15.9800
   14.2700
   15.8100
   15.5200
   21.6500
   14.5400
   19.0000
   16.1400
   15.9500
   15.9000
   17.9100
   14.3000
   19.9400
   13.8600
   14.4400
   16.0000
   18.3200
   10.0000
   18.1700
   16.8700
   13.9800
   20.2900
   10.2500
   18.0000
   21.6700
   20.9600
   16.1600
   15.7800
   15.9600
   16.1700
    9.1940
   13.8100
   19.3400];

%% Clear temporary variables
clear opts


%% sorting/naming data
Monkey=table2array(DurationNeuroS1(:,1);
Locs=table2array(DurationNeuroS1(:,4));
Thresholds=table2array(DurationNeuroS1(:,9));
cfs=table2array(DurationNeuroS1(:,6));
durs=table2array(DurationNeuroS1(:,10));
Nthresholds=table2array(DurationNeuroS1(:,12));
N2threshold=table2array(DurationNeuroS1(:,13));
Nslope=table2array(DurationNeuroS1(:,14));
N2slope=table2array(DurationNeuroS1(:,15));

%% gets thresholds, fits with exp, gets time constant and const. of Prop.
unitct=1; % keeping track of what unit we're on

for threshct = 1:7:length(Thresholds)

   
    
    % variables that will be fed to fitnlm
    durscheck(threshct:threshct+6,1)=durs(threshct:threshct+6,1);
    y=Thresholds(threshct:threshct+6,1);
    
    if isnan(y(1,1))
        continue
    end
    
    x=[200,100,50,25,12.5,6.5,3.25];
    
    % fit the data
    tauguess=1;
    rangeguess=0.01;
    mdl=taufxn_v2(x,y,tauguess,rangeguess); % function that fits with exponential
    % consider adding in if statement here to fix fits
    tau(unitct,2)=mdl(1,2); % extract tau
    tau(unitct,4)=mdl(1,1);%const of prop
    tau(unitct,1)=cfs(threshct,1)/1000; %putting cfs in kilohertz
    tau(unitct,3)=Locs(threshct,1); %keeping track of brain region
    
    
    unitct=unitct+1;
end


%% gets time constants for noise. i know i know, shouldn't copy/paste like this
unitct=1; %reset the unit count
% the noise thresholds
for noisethreshct = 1:7:length(Nthresholds)

   
    
    % variables that will be fed to fitnlm
%     durscheck(noisethreshct:noisethreshct+6,1)=durs(noisethreshct:noisethreshct+6,1);
    y=Nthresholds(noisethreshct:noisethreshct+6,1);
    
    if isnan(y(1,1))
        continue
    end
    
    x=[200,100,50,25,12.5,6.5,3.25];
    
    % fit the data
    tauguess=1;
    rangeguess=0.01;
    mdl=taufxn_v2(x,y,tauguess,rangeguess); % function that fits with exponential
    % consider adding in if statement here to fix fits
    ntau(unitct,2)=mdl(1,2); % extract tau
    ntau(unitct,4)=mdl(1,1);%const of prop
    ntau(unitct,1)=cfs(noisethreshct,1)/1000; %putting cfs in kilohertz
    ntau(unitct,3)=Locs(noisethreshct,1); %keeping track of brain region
    
    
    unitct=unitct+1;
end

unitct=1;
% 44 dB noise
for noise44threshct = 1:7:length(N2threshold)

   
    
    % variables that will be fed to fitnlm
%     durscheck(noisethreshct:noisethreshct+6,1)=durs(noisethreshct:noisethreshct+6,1);
    y=N2threshold(noise44threshct:noise44threshct+6,1);
    
    if isnan(y(1,1))
        continue
    end
    
    x=[200,100,50,25,12.5,6.5,3.25];
    
    % fit the data
    tauguess=1;
    rangeguess=0.01;
    mdl=taufxn_v2(x,y,tauguess,rangeguess); % function that fits with exponential
    % consider adding in if statement here to fix fits
    n44tau(unitct,2)=mdl(1,2); % extract tau
    n44tau(unitct,4)=mdl(1,1);%const of prop
    n44tau(unitct,1)=cfs(noise44threshct,1)/1000; %putting cfs in kilohertz
    n44tau(unitct,3)=Locs(noise44threshct,1); %keeping track of brain region
    
    
    unitct=unitct+1;
end

%% plotting time constants

% quiet
% subplot(1,3,1)
% gscatter(tau(:,1),tau(:,2),tau(:,3))
% title('quiet')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% axis square
% legend('CN','IC')
% ylabel('Time Constant (ms)')
% xlabel('Tone Frequency (kHz)')
% 
% % 34 db noise
% subplot(1,3,2)
% gscatter(ntau(:,1),ntau(:,2),ntau(:,3))
% title('34 dB noise')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% axis square
% legend('CN','IC')
% xlabel('Tone Frequency (kHz)')
% 
% % 44 db noise
% subplot(1,3,3)
% gscatter(n44tau(:,1),n44tau(:,2),n44tau(:,3))
% title('44 dB noise')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% axis square
% legend('CN','IC')
% xlabel('Tone Frequency (kHz)')
% 
% %setting up for the histograms below
% ICidx=find(tau(:,3)==2);
% CNidx=find(tau(:,3)==1);
% CNidxN=find(n44tau(:,3)==1);
% ICidxN=find(n44tau(:,3)==2);
% edges=0:5:150;
% 
% 
% %histograms, which show how IC and CN are diff

figure
histogram(tau(ICidx,2),edges)
hold on
histogram(tau(CNidx,2),edges)
hold on
%histogram(tau_beh,edges)
% legend('IC','CN')
title('Time Const. in Q')
ylabel('# Units')
xlabel('Time Constant')
% 


% figure
% histogram(tau(CNidx,2),edges)
% hold on
% histogram(n44tau(CNidxN,2),edges)
% legend('quiet','noise')
% title('CN time constants')
% ylabel('# Units')
% xlabel('Time Constant')
% 
% figure
% histogram(tau(ICidx,2),edges)
% hold on
% histogram(n44tau(ICidxN,2),edges)
% legend('quiet','noise')
% title('IC time constants')

figure
cdfplot(tau(CNidx,2))
hold on
cdfplot(n44tau(CNidxN,2))
hold on
cdfplot(tauN_beh)
legend('Quiet','Noise','Behavior')
title('CN time constants')

figure
cdfplot(tau(ICidx,2))
hold on
cdfplot(n44tau(ICidxN,2))
hold on
cdfplot(tauN_beh)
legend('Quiet','Noise','Behavior')
title('IC time constants')
% 
% figure
% histogram(n44tau(ICidxN,2),edges)
% hold on
% histogram(n44tau(CNidxN,2),edges)
% hold on
% histogram(tauN_beh,edges)
% legend('IC','CN','Behavior')
% title('Time constants in Noise')
% ylabel('# Units')
% xlabel('Time Constant')

% figure
% histogram(n44tau(ICidxN,2),edges)
% hold on
% histogram(n44tau(CNidxN,2),edges)
% hold on
% histogram(tauN_beh,edges)
% title('Noise')
% ylabel('# Units')
% xlabel('Time Constant (ms)')


% 
% 
% %% plotting constant or proportionality
% figure
% 
% % const of prop in quiet
% subplot(1,3,1)
% gscatter(tau(:,1),tau(:,4),tau(:,3))
% title('quiet')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% axis square
% legend('CN','IC')
% ylabel('Const. of Prop')
% xlabel('Tone Frequency (kHz)')
% 
% % const of prop in 34
% subplot(1,3,2)
% gscatter(ntau(:,1),ntau(:,4),ntau(:,3))
% title('34 dB Noise')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% axis square
% legend('CN','IC')
% ylabel('Const. of Prop')
% xlabel('Tone Frequency (kHz)')
% 
% % const of prop in 44
% subplot(1,3,3)
% gscatter(n44tau(:,1),n44tau(:,4),n44tau(:,3))
% title('44 dB Noise')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% axis square
% legend('CN','IC')
% ylabel('Const. of Prop')
% xlabel('Tone Frequency (kHz)')



%% stats - komogorov smirnov to compare IC, CN and Behavior in Q and N

% [h,p] = kstest2(tau(ICidx,2),tau(CNidx,2)) % cn and ic are not diff in quiet (p = 0.29)
% [hnoise,p] = kstest2(n44tau(ICidxN,2),n44tau(CNidxN,2)) % ic and cn are diff in noise (p = 0.022)
[hh,pval] = kstest2(tau(ICidx,2),n44tau(ICidxN,2)) % ic is faster in quiet than in noise(p = 0.017)
% [h_beh,p_beh] = kstest2(tau_beh,tau(ICidx,2))
% [h_behN,p_behN] = kstest2(tauN_beh,tau(ICidxN,2))
% [h_behCNq,pbehCNq] = kstest2(tau_beh,tau(CNidx,2))
% [h_behCNn,pbehCNn] = kstest2(tauN_beh,tau(CNidxN,2)) %all the neuro is diff from behavior

% also want to just compare thresholds and slopes


%% mixed effects model
t=table(Locs,Thresholds,N2threshold,slopes,N2slope,cfs,durs,Monkey);
lmem_threshQ=fitlme(t,'Thresholds~durs*Locs+cfs+(1|Monkey)')
lmem_slopeQ=fitlme(t,'slopes~Locs+durs*cfs+(1|Monkey)')
lmem_threshN=fitlme(t,'N2threshold~durs*Locs+cfs+(1|Monkey)')
lmem_slopesN=fitlme(t,'N2slope~durs*Locs+cfs+(1|Monkey)')
















