%% Modulation Transfer function "within cell" pooling
% MTF analysis/plotting written by Sam Hauser (~2015), pooling and ROC analysis written by
% Chase Mackey (2021)
%
% this is the version that loops through multiple blocks 
% and pools responses pre-analysis like johnson et al. (2012, j neurophys)
% "within cell" method

% this script is importantly different than samsmoodlongpoolerfxn and
% samsmodlongpoolerlooper.m, which (together) sample single-trials and do
% the johnson et al. (2012, j neurophys) "across cell" method
%
clear all; 
close all;

%% set variables

%respNonPref = onefreq(1,19); % setting the comparison stimulus
stdidx=20; % comparison stim for rate
bmf=24; %bmf for rate
stdvs=5; %comparison stim for VS 
bmfvs=10; %bmf for vs

flipr=0; %flip the sign if you're using downward sloping part of MTF
flipvs=0; %flip sign for vector strength

% duration
int = 2; %interval for stepping through binned spikes (2 for 500 ms, 4, for 250 ms, etc.)
dur = 1/int; %bin width in seconds (0.5, 0.25, 0.1, 0.05)

% make this 1 if you want to do ROC analysis
discanalysis=1;

%% import all the data

paths = {
        
        '/Volumes/ramlab/Charlie data/Neurobehavior/IC/ex180308/tank180308c/OurData-12';
        %'/Volumes/ramlab/Delta data/Neurobehavior/IC/ex150727/tank150727d/OurData-14';
        %'/Volumes/ramlab/Delta data/Neurobehavior/IC/ex150917/tank150917d/OurData-4';
        %'/Volumes/ramlab/Delta data/Neurobehavior/IC/ex151012/tank151012d/OurData-19';
        
        };
    
    

check = 0;
SortCode=0; % sort code was usually zero
for ct = 1:1:size(paths)
    
    %import
    data = TDTbin2mat(paths{ct,:}); 
    allspikes = [];
    allspikes = data.snips.eNeu.ts;
    
    %sort spikes using sort code
    n=1;
    for i = 1:length(data.snips.eNeu.sortcode)
        if data.snips.eNeu.sortcode(i,1) == SortCode 
            sortedspikestmp(n,1) = allspikes(i,1);
            n = n+1;
        end
    end
    
    %make an array with spikes from multiple blocks of data
    if check==0
        sortedspikes = sortedspikestmp;
        check = check+1;
    else
       sortedspikes=cat(1,sortedspikes,sortedspikestmp); 
    end
end


%% SORT SPIKES

% To get just the ones from the newly sorted file. If you don't want the
% sorted one, you can change the sortcode number to 0 or you can change all
% of the places where it says sortedspikes back to allspikes below this. 


% n=1;
% sortedspikes = []; 
% for i = 1:length(data.snips.eNeu.sortcode)
%     if data.snips.eNeu.sortcode(i,1) == SortCode 
%         sortedspikes(n,1) = allspikes(i,1);
%         n = n+1;
%     end
% end


%% old way of binning %%
% bins = 0.25:0.5:116.25; %This centers the bins over the right one second interval
% figure; 
% hist(sortedspikes,bins); %Let's look at it, then save that list
% binnedspc = hist(sortedspikes,bins);
%
%% new way of binning where duration is manipulated (uses bin edges not center)
% figure
edges=0:dur:116.5; % left bin edges
figure
binnedspc=histogram(sortedspikes,edges)

% now try to manipulate bin width/duration of window where we count spikes
% figure
% histogram(sortedspikes,'BinWidth',0.5,'BinEdges',edges)

% assign spikes to this variable
binnedspc=binnedspc.Values;

%% COUNT SPIKES
stimspc = zeros(116,1); %This is just the bins that have the stimulus in it. 
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


%First get timestamps lined up
spikesortbysec = [];
col = 1;
for i = 1:int:length(binnedspc)
    if i > 1
        sofar = sum(binnedspc(1, 1:i-1));
        spikesortbysec(1:binnedspc(1, i), col) = sortedspikes(sofar+1:sofar+binnedspc(1,i), 1) ;
        col = col + 1;
    else 
        spikesortbysec(1:binnedspc(1, i), col) = sortedspikes(1:binnedspc(1,i), 1);
        col = col + 1;
    end
end

% collapse frequency
% spikebyfreq = []; 
% n = 1;
% for i = 1:116
%     
%     spikebyfreq(1:sum(stimspc(i:(i+1), 1)), n) = vertcat(spikesortbysec(1:stimspc(i, 1), i), spikesortbysec(1:stimspc((i+1), 1), (i+1)));
%     n = n + 1; 
% end
%%        

% then run this to find the VS at each frequency
r = zeros(116,1); 
radpersec = zeros(1, length(freq));


for Z = 1:length(freq)
    radpersec(1, Z) = (2*pi)*freq(1,Z);
    nZ=stimspc(Z,1);
    reps=spikesortbysec(1:nZ,Z);
    phaseangleZ=radpersec(1, Z) .* reps;
    xZ=cos(phaseangleZ);
    yZ=sin(phaseangleZ);
    r(Z,1)=sqrt((sum(xZ))^2+(sum(yZ))^2)/nZ;
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
            whichspikeVS(VSsort,column) = VSsort2;
        end
    end
end

OUT = [ onefreq' avestimspc' done];

%% Neuronal discrimination analysis
if discanalysis ==1
    %% ROC analysis on spike count

    [pc,FA,HR]=mfROC(levelspikes,stdidx,bmf,0); % pass data to ROC function
    pc(:,2) = onefreq(1,stdidx:bmf)'; % get the mod freqs


    %% ROC analysis on vector strength 

     [PCvs,FAvs,HRvs]=mfROC(levelspikesVS(1:end-1,:),stdvs,bmfvs,1); % pass data to ROC function
     PCvs(:,2) = onefreq(1,stdvs:bmfvs)'; % get the mod freqs

    %% plot neurometric functions with MTF data
    
    % fit with weibull to extract thresh and slope
   [ratethresh,ratethreshdelta,sloper]=createFitMF(pc(:,2),pc(:,1),pc(1,2),0,flipr);
   
   [VSthresh,VSthreshdelta,slopevs]=createFitMF(PCvs(:,2),PCvs(:,1),PCvs(1,2),0,flipvs); 
   
   % printing these out
   slopes=[sloper,slopevs];
   DiscThresholds = [ratethresh,VSthresh,ratethreshdelta,VSthreshdelta, sloper,slopevs]
   
   if flipvs == 1
       pcvsflip = PCvs;
       pcvsflip(:,1) = 1-pcvsflip(:,1);
       linthreshVS=ThreshCalc(pcvsflip')-pcvsflip(1,2)
   end
   
   if flipr == 1
       pcrflip = pc;
       pcrflip(:,1) = 1-pcrflip(:,1);
       linthreshRate=ThreshCalc(pcrflip')-pcrflip(1,2)
   end
   
   % PLOTTING EVERYTHING TOGETHER
    subplot(2,2,1)
    semilogx(onefreq,avestimspc, onefreq, SSave, 'r*'); 
    title('Spike Count');
    xlabel('Modulation Frequency'); 
    ylabel('Firing Rate (Hz)')
    subplot(2,2,2)
    semilogx(onefreq, done);
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
