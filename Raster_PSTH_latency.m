%% Raster Peristimulus Time Histogram (PSTH) Example
%
% <html>
% Import snippet and epoc data into Matlab using TDTbin2mat <br>
% Generate peristimulus raster and histogram plots over all trials <br>
% Good for stim-response experiments, such as optogenetic or electrical
% stimulation
% </html>


close all; clear all; clc;

Monkey='Delta';
Date='121017d';
Num='9';
Loc='IC';
Num=num2str(Num);
Date1=Date(1:6);
MyTank = ['Y:\' Monkey ' data\NeuroBehavior\' Loc '\ex' Date1 '\tank' Date];
MyBlock = ['~OurData-' Num];
BLOCKPATH = ('/Volumes/ramlab/Charlie data/Neurobehavior/IC/ex180307/tank180307c/OurData-27');
threshV=2.50000e-04;
%%
% Set up the variables for the data you want to extract. We will extract
% channel 1 from the eNe1 snippet data store, created by the PCA Sorting
% gizmo, and use our PulseGen epoc event ('PC0/') as our stimulus onset.
REF_EPOC = 'Freq';
SNIP_STORE = 'eNeu';
SORTID = 'TankSort';
CHANNEL = 1;
SORTCODE = 0; % set to 0 to use all sorts
TRANGE = [-0.25, 0.45]; % window size [start time relative to epoc onset, window duration]

%%
% Now read the specified data from our block into a Matlab structure. The
% 'NODATA' flag means that we are only intereseted in the snippet
% timestamps, not the actual snippet waveforms in this example.
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'snips', 'scalars'}, 'SORTNAME', SORTID, 'CHANNEL', CHANNEL, 'NODATA', 1);


%% Use TDTfilter to extract data around our epoc event
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event.
raster_data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE);

%%
% Adding the 'TIMEREF' flag makes all of the timestamps relative to the
% epoc event, which is ideal for generating histograms.
hist_data = TDTfilter(data, REF_EPOC, 'TIME', TRANGE, 'TIMEREF', 1);

%% 
%data is now in Matlab. The rest of the code is a simple plotting example.  
%First, we'll find matching timestamps for our 
% selected sort code (unit).
TS = raster_data.snips.(SNIP_STORE).ts;
if SORTCODE ~= 0
    i = find(raster_data.snips.(SNIP_STORE).sortcode == SORTCODE);
    TS = TS(i);
end
if isempty(TS)
    error('no matching timestamps found')
end

num_trials = size(raster_data.time_ranges, 2);

%% Make the histogram plot
figure('Position',[100, 100, 500, 800]);

hist_TS = hist_data.snips.(SNIP_STORE).ts;
subplot(2,1,1);
NBINS = floor(numel(hist_TS)/10);
hist(hist_TS, NBINS);
N = hist(hist_TS, NBINS); hold on;

axis tight; axis square;
set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
ylabel('Count','FontSize',16);
title({'Peristimulus', sprintf('Channel %d, n = %d trials', CHANNEL, num_trials)});

% Draw a vertical line at t=0.
line([0 0], [0, max(N)], 'Color','r', 'LineStyle','-', 'LineWidth', 3);

% Creating the Raster Plot
% For the raster plot, make a cell array of timestamps for each trial.
all_TS = cell(num_trials, 1);
all_Y = cell(num_trials, 1);
for trial = 1:num_trials
    trial_on = raster_data.time_ranges(1, trial);
    trial_off = raster_data.time_ranges(2, trial);
    trial_TS = TS(TS >= trial_on & TS < trial_off);
    all_TS{trial} = trial_TS - trial_on + TRANGE(1);
    all_Y{trial} = trial * ones(numel(trial_TS), 1);
end
all_X = cat(1, all_TS{:});
all_Y = cat(1, all_Y{:});

% Make the raster plot.
% subplot(2,1,2);
% plot(all_X, all_Y, '.', 'MarkerEdgeColor','k', 'MarkerSize',15); hold on;
% axis tight; axis square;
% set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
% xlabel('Trial Window, s','FontSize',16);
% ylabel('Trial Number','FontSize',16);
% title({'Raster', sprintf('Channel %d, n = %d trials', CHANNEL, num_trials)});
% 
% % Draw a vertical line at t=0.
% line([0 0], [0, trial+1], 'Color','r', 'LineStyle','-', 'LineWidth', 3);


%% Extracting spike times at a single lvl for CDF/latency calc.

% % get the levels used 
[freq levl nlvl nLvls Responses TT figure1] = DataExtraction(MyTank,MyBlock);
[lvls,alpha]=TLevlFindNOrder(levl);
close figure 1

% round
cdflvl=dBtoSPL(threshV,1); % convert it to dB
dec=1e6; % round to 5th decimal
levl_cdf = sym(levl(1,:),'d'); % convert to sym because we want to round
levl_cdf = round(levl_cdf*dec);
levl_cdf=double(levl_cdf/dec); % convert back to double so it's not useless


% % find every trial where lvl == threshV and extract those trials
cdfSpikeTimes={};
for tr_ct=1:length(levl_cdf) % y was previously defined as the num of trials
    if levl_cdf(1,tr_ct)==threshV
       cdfSpikeTimes(tr_ct,1)=all_TS(tr_ct,1);
    end
end
% now cdfSpikeTimes has the spike times from the right trials

%removing empties
cdfSpikeTimes(cellfun('isempty',cdfSpikeTimes))=[];

%% compute baseline activity
   
%collapse across trials and sort in ascending order
 sortedtimes=cell2mat(cdfSpikeTimes);
 sortedtimes=sort(sortedtimes);
 
% number of pre stim spikes
 negs=1;
 for sortct=1:length(sortedtimes)
     if sortedtimes(sortct,1)<=0
         negs=negs+1;
     end
 end
 spont=negs/length(cdfSpikeTimes); %avg spikes in 250ms before tone
 spontR=(spont*4)/1000;
 

 
 
%% Cumulative spike fxn
% consider spikes from a certain time range
pre=-0.25;
post=0.2;

% making a csf for every trial
spkcum=[];
cdfbin = ((pre*1000):(post*1000))/1000;
for ii = 1:length(cdfSpikeTimes)
    for jj = 1:length(cdfbin)
        spkcum(ii,jj) = sum((cdfSpikeTimes{ii}<cdfbin(jj)));
    end
%     spontR=spontR+spontR;
end

% cumulative spont activity estimate
cumspont=zeros(length(cdfSpikeTimes),1);
for kk=1:length(cdfSpikeTimes)
    cumspont(kk,1)=cumspont(kk,1)+SpontR;
end

% subtracting spont from cumulative spike functions (Rowland et al. 2007)
for ff=1:length(spkcum)
    spkcum(ff,1)=spkcum(ff,1)-cumspont(ff,1);
end
%spkcum(:,:)=spkcum(:,:)-spontR;


% avg across bins to get one csf
avgspkcum=[];
for avgct=1:length(spkcum)
    avgspkcum(avgct,1)=mean(spkcum(:,avgct));
end




% trying to plot polyfit
% xlinfit(1,1)=min(cdfbin);
% xlinfit(2,1)=max(cdfbin);
% ylinfit(1,1)=min(avgspkcum);
% ylinfit(2,1)=max(avgspkcum);
% 
% p=polyfit(xlinfit,ylinfit,2);
% 
% r = p(1) .* cdfbin(:,1) + p(2);
% hold on
% plot (cdfbin,avgspkcum,r)

%% plot cumulative spike function
plot(cdfbin,avgspkcum)
hold on
figure
cdfplot(sortedtimes)