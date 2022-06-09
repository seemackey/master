%% Modulation Transfer function analysis
% MTF analysis/plotting written by Sam Hauser (~2015), ROC analysis written by
% Chase Mackey (2021)

clear all; 
close all;

%% import 
data = TDTbin2mat('/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151218/tank151218c/OurData-13'); % trying the import without active x 
%data = TDTbin2mat('/Volumes/ramlab/Delta data/Neurobehavior/IC/ex151120/tank151120d/OurData-7'); % trying the import without active x 

%% set variables
SortCode=0; % sort code was usually zero
%respNonPref = onefreq(1,19); % setting the comparison stimulus
stdidx=15; % comparison stim for rate
bmf=28; %bmf for rate
stdvs=22; %comparison stim for VS 
bmfvs=28; %bmf for vs

flipr=0; %flip the sign if you're using downward sloping part of MTF
flipvs=1; %flip sign for vector strength

dur = 0.5; %bin width in seconds (0.5, 0.25, 0.1, 0.05)
int = 2; %interval for stepping through binned spikes (2 for 500 ms, 4, for 250 ms, etc.)

% make this 1 if you want to do ROC analysis
discanalysis=1;

%% SORT SPIKES

% To get just the ones from the newly sorted file. If you don't want the
% sorted one, you can change the sortcode number to 0 or you can change all
% of the places where it says sortedspikes back to allspikes below this. 

allspikes = data.snips.eNeu.ts;
n=1;
sortedspikes = []; 
for i = 1:length(data.snips.eNeu.sortcode)
    if data.snips.eNeu.sortcode(i,1) == SortCode % 3 because it's sortcode is under the number 3
        sortedspikes(n,1) = allspikes(i,1);
        n = n+1;
    end
end


%% changing length of analysis window based on mod. freq
% modfreq=[64	25.39841683	20.1587368	812.7493386	1024	2	203.1873347	101.5936673	256	12.69920842	40.3174736	6.349604208	3.174802104	2.5198421	5.0396842	10.0793684	16	322.5397888	4	32	50.79683366	80.63494719	8	0	128	406.3746693	161.2698944	512	645.0795775	161.2698944	6.349604208	32	0	512	64	5.0396842	40.3174736	1024	203.1873347	20.1587368	80.63494719	8	812.7493386	2	50.79683366	16	406.3746693	10.0793684	256	128	101.5936673	645.0795775	4	25.39841683	12.69920842	322.5397888	3.174802104	2.5198421	50.79683366	80.63494719	12.69920842	40.3174736	32	161.2698944	10.0793684	512	645.0795775	16	322.5397888	812.7493386	20.1587368	1024	5.0396842	6.349604208	128	256	8	0	4	2.5198421	25.39841683	406.3746693	64	3.174802104	2	101.5936673	203.1873347	50.79683366	128	5.0396842	256	64	10.0793684	645.0795775	20.1587368	80.63494719	512	40.3174736	12.69920842	16	2.5198421	6.349604208	2	101.5936673	8	406.3746693	0	25.39841683	161.2698944	322.5397888	203.1873347	32	4	812.7493386	1024	3.174802104];
% modfreq=roundn(modfreq,-7);
% 
% %calculates time-windows (modfreq(trim(3,:)) so no resps to incomp. AM are used
% [modfreq_trim] = trimmer(modfreq,dur*1000);
% 
% modfreq_trim(3,:)=modfreq_trim(3,:)/1000;
% 
% % counting spikes
% for ts_ct = 1:1:length(modfreq_trim)
%    stimspc(ts_ct,1) = length(find(sortedspikes(:,1)<=((ts_ct-1)+modfreq_trim(3,ts_ct))&(sortedspikes(:,1))>ts_ct-1));
% end


%% old way of binning %%
% bins = 0.25:0.5:116.25; %This centers the bins over the right one second interval
% figure; 
% hist(sortedspikes,bins); %Let's look at it, then save that list
% binnedspc = hist(sortedspikes,bins);
%
%% new way of binning where duration is manipulated (uses bin edges not center)
% figure
% 
% dur_alt=modfreq_trim(3,:)/1000;

edges=0:dur:116.5; % left bin edges
figure
binnedspc=histogram(sortedspikes,edges);


% now try to manipulate bin width/duration of window where we count spikes
% figure
% histogram(sortedspikes,'BinWidth',0.5,'BinEdges',edges)

% assign spikes to this variable
binnedspc=binnedspc.Values;



%% COUNT SPIKES
stimspc = zeros(116,1); %will be the bins that have the stimulus in it. 
n = 0;
for i = 1:int:length(binnedspc)
    n = n+1;
    stimspc(n,1) = binnedspc(1,i);
end

avestimspc = zeros(29,1);

n = 0; 
onefreq = 2.^[1:(1/3):10];
onefreq(1, 29) = 0; 
onefreq = roundn(onefreq,-7); 
modfreq=[64	25.39841683	20.1587368	812.7493386	1024	2	203.1873347	101.5936673	256	12.69920842	40.3174736	6.349604208	3.174802104	2.5198421	5.0396842	10.0793684	16	322.5397888	4	32	50.79683366	80.63494719	8	0	128	406.3746693	161.2698944	512	645.0795775	161.2698944	6.349604208	32	0	512	64	5.0396842	40.3174736	1024	203.1873347	20.1587368	80.63494719	8	812.7493386	2	50.79683366	16	406.3746693	10.0793684	256	128	101.5936673	645.0795775	4	25.39841683	12.69920842	322.5397888	3.174802104	2.5198421	50.79683366	80.63494719	12.69920842	40.3174736	32	161.2698944	10.0793684	512	645.0795775	16	322.5397888	812.7493386	20.1587368	1024	5.0396842	6.349604208	128	256	8	0	4	2.5198421	25.39841683	406.3746693	64	3.174802104	2	101.5936673	203.1873347	50.79683366	128	5.0396842	256	64	10.0793684	645.0795775	20.1587368	80.63494719	512	40.3174736	12.69920842	16	2.5198421	6.349604208	2	101.5936673	8	406.3746693	0	25.39841683	161.2698944	322.5397888	203.1873347	32	4	812.7493386	1024	3.174802104];
modfreq=roundn(modfreq,-7);
freq=modfreq;
whichspike = zeros(29,5); 
whichspike(:,5) = onefreq';

levelspikes = zeros(29,4); 
for i = 1:29
    col = 0; 
    for n = 1:116
        if modfreq(1, n) == onefreq(1,i)
            col = col + 1; 
            levelspikes(i,col) = stimspc(n, 1);
            whichspike(i,col) = n;
        end
    end
end
%%   
avestimspc = nanmean(levelspikes'); 
avestimspcstd = std(levelspikes'); 

y = avestimspc(1, 29); 

for i = 1:29 
    SSave(1, 1:i) = y;
end


%Just looking at the rMTF
if discanalysis ==0
    figure; 
    semilogx(onefreq,avestimspc, '*r', onefreq, SSave, '-');
    figure; 
    errorbar(avestimspc, avestimspcstd, 'r*'); 

    title('Spike Count');
    xlabel('Modulation Frequency'); 
    ylabel('# of Spikes per second')
end

%% Vector Strength! 

%Vector strength timing analysis script.
%for the current application, this will take the modulated noise
%signal and calculate the vector strength of the spikes from physiology
%according to the methods in Goldberg et al 1969.
%theta=phase angle
%x=cos(theta)
%y=sin(theta)
%direction is a measure of mean phase relation between the stimulus and
%discharge given by: phi=arctan(sum(y)/sum(x))+k*(pi), where k is zero or
%one depending on signs of sum(x) and sum(y)
%vector strength, r, is the length of the mean vector.
%r=squareroot(sum(x)^2+sum(y)^2)/n
%n=number of vectors

% March 2022 - added phase projected vector strength


%First get timestamps lined up
spikesortbysec = [];
col = 1;
for i = 1:2:length(binnedspc)
    if i > 1
        sofar = sum(binnedspc(1, 1:i-1));
        spikesortbysec(1:binnedspc(1, i), col) = sortedspikes(sofar+1:sofar+binnedspc(1,i), 1) ;
        col = col + 1;
    else 
        spikesortbysec(1:binnedspc(1, i), col) = sortedspikes(1:binnedspc(1,i), 1);
        col = col + 1;
    end
end

%collapse frequency
spikebyfreq = []; 
n = 1;
for i = 1:116
    
    spikebyfreq(1:sum(stimspc(i:(i+1), 1)), n) = vertcat(spikesortbysec(1:stimspc(i, 1), i), spikesortbysec(1:stimspc((i+1), 1), (i+1)));
    n = n + 1; 
end
%%        

% then run this to find the VS at each frequency
r = zeros(116,1); 
radpersec = zeros(1, length(freq));
angles={length(freq),2};

for Z = 1:length(freq)
    radpersec(1, Z) = (2*pi)*freq(1,Z);
    nZ=stimspc(Z,1);
    reps=spikesortbysec(1:nZ,Z);
    phaseangleZ=radpersec(1, Z) .* reps;
    angles{Z,1}=phaseangleZ;
    xZ=cos(phaseangleZ);
    yZ=sin(phaseangleZ);
    phi(Z,1)=atan2(sum(yZ),sum(xZ)); % used for VSpp
    r(Z,1)=sqrt((sum(xZ))^2+(sum(yZ))^2)/nZ;    % classical VS
end



avvs=zeros(29,1);
done= [];
width=size(levelspikes);
for i = 1:29
    for n = 1:116
        if onefreq(1,i) == modfreq(1, n)
            avvs(i,1) = avvs(i,1) + r(n, 1); 
        end
    end
    done(i,1) = avvs(i,1)/4; % need to make this so the 4 isn't hard coded
end

% just looking at the vsMTF
if discanalysis==0
    figure; 
    semilogx(onefreq, done);
    title('Vector Strength'); 
    xlabel('Modulation Frequency');
    ylabel('Vector Strength'); 
end

%% sort VS by mod freq

levelspikesVS = zeros(29,4); 
for VSsort = 1:length(onefreq)
    column = 0; 
    for VSsort2 = 1:length(freq)
        if modfreq(1, VSsort2) == onefreq(1,VSsort)
            column = column + 1; 
            levelspikesVS(VSsort,column) = r(VSsort2, 1);
            phi_sort(VSsort,column) = phi(VSsort2,1);
            angles_sort{VSsort,column} = angles{VSsort2,1};
            whichspikeVS(VSsort,column) = VSsort2;
        end
    end
end

% phase projected vector strength (Yin et al. 2011 j neurophys) even though
% it's stupid
for trial = 1:1:length(phi_sort(1,:))
    for mf = 1:1:length(phi_sort(:,1))
        angles_tmp = vertcat(angles_sort{mf,1:4});
        phi_c(mf,1) = atan2(sum(sin(angles_tmp)),sum(cos(angles_tmp)));
        VSpp(mf,trial) = levelspikesVS(mf,trial)*cos(phi_sort(mf,trial)-phi_c(mf,1));
    end
end

% atan2(sum(((length(reps)).*(sin(angles{Z,1})))),sum(((length(reps)).*(cos(angles{Z,1})))))

OUT = [ onefreq' avestimspc' done];

%% Neuronal discrimination analysis
if discanalysis ==1
    %% ROC analysis on spike count

    [pc,FA,HR]=mfROC(levelspikes,stdidx,bmf,0); % pass data to ROC function
    pc(:,2) = onefreq(1,stdidx:bmf)'; % get the mod freqs


    %% ROC analysis on vector strength 

     %[PCvs,FAvs,HRvs]=mfROC(levelspikesVS(1:end-1,:),stdvs,bmfvs,1); % pass data to ROC function
     [PCvs,FAvs,HRvs]=mfROC(VSpp(1:end-1,:),stdvs,bmfvs,1); % pass data to ROC function
     PCvs(:,2) = onefreq(1,stdvs:bmfvs)'; % get the mod freqs

    %% plot neurometric functions with MTF data
    
    % fit with weibull to extract thresh and slope
   %[ratethresh,ratethreshdelta,sloper]=createFitMF(pc(:,2),pc(:,1),pc(1,2),0,flipr);
   
   % flip data if we need to 
   if flipr==0
   [linthresh] = ThreshCalc(pc');
   elseif flipr==1
       pc(:,1)=1-pc(:,1);
   [linthresh] = ThreshCalc((pc'));
   end
   
   % linear interp threshold
   linthresh=linthresh-pc(1,2)
   
   % flip data if we need to 
   if flipvs==0
       [linthreshVS] = ThreshCalc(PCvs');
   elseif flipvs==1
       PCvs(:,1)=1-PCvs(:,1);
   [linthreshVS] = ThreshCalc((PCvs'));
   end
   
   % linear interp threshold
   linthreshVS
   linthreshVS=linthreshVS-PCvs(1,2)
   
   %[VSthresh,VSthreshdelta,slopevs]=createFitMF(PCvs(:,2),PCvs(:,1),PCvs(1,2),0,flipvs); 
   
   % printing these out
   %slopes=[sloper,slopevs];
   %DiscThresholds = [ratethresh,VSthresh,ratethreshdelta,VSthreshdelta, sloper,slopevs];
   
   % PLOTTING EVERYTHING TOGETHER
   figure
    subplot(2,2,1)
    semilogx(onefreq,avestimspc, onefreq, SSave, 'r*'); 
    title('Spike Count');
    xlabel('Modulation Frequency'); 
    ylabel('Firing Rate (Hz)')
    subplot(2,2,2)
    semilogx(onefreq, done);
    hold on
    semilogx(onefreq,mean(VSpp,2));
    legend('VS','VSpp')
    title('Vector Strength'); 
    xlabel('Modulation Frequency');
    ylabel('Vector Strength'); 
    subplot(2,2,3)
    plot(pc(:,2),pc(:,1),'-o') % plot probability correct in a 2AFC
    title('Rate-based neurometric function'); 
    xlabel('Modulation Frequency');
    ylabel('Neurometric P(C)');
    subplot(2,2,4)
    plot(PCvs(:,2),PCvs(:,1),'-o') % plot probability correct in a 2AFC
    title('VS-based neurometric function');
    xlabel('Modulation Frequency');
    ylabel('Neurometric P(C)');
    
    % saving these for manual curve fitting
     x=[PCvs(:,2)]; % vector strength
     y=[PCvs(:,1)];
    
    xr=[pc(:,2)]; % rate
    yr=[pc(:,1)];
end
