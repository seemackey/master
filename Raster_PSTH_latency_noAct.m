%% Raster Peristimulus Time Histogram
% with cumulative spike fxn for latency calc.
% Chase Mackey, 2019


close all; clear all; clc;

%% ENTER PATH INFO HERE
Monkey='Charlie';
Date='130114c';
Num='6';
Loc='IC';
Num=num2str(Num);
Date1=Date(1:6);

MyTank = ['Y:\' Monkey ' data\NeuroBehavior\' Loc '\ex' Date1 '\tank' Date];
MyBlock = ['~OurData-' Num];
%BLOCKPATH = ['Y:\' Monkey ' data\NeuroBehavior\' Loc '\ex' Date1 '\tank' Date '\OurData-' Num];
BLOCKPATH='/Volumes/ramlab/Charlie data/Neurobehavior/IC/ex130114/tank130114c/OurData-6'
%cf=18000;
% ENTER LVL FOR LATENCY CALC.
% threshV=0.003530000103638;
% %threshV=[2.499999936844688e-02];
% %threshV=[0.0249999994412065];
% threshdB=dBtoSPL(threshV,1)
% dec=1e5; % round to nth decimal later in script
%%%%%%%%
%%
% Set up the variables for the data you want to extract. We will extract
% channel 1 from the eNeu snippet data store, sometimes snipstore = Snip,
% sometimes it's eNeu
REF_EPOC = 'Freq';
SNIP_STORE = 'eNeu'; %this variable is called eNeu after ~2011 or so
SORTID = 'TankSort';
CHANNEL = 1;
SORTCODE = 0; % set to 0 to use all sorts
TRANGE = [-0.25, 0.6]; % window size [start time relative to epoc onset, window duration]

%%
% Now read the specified data from our block into a Matlab structure. The
% 'NODATA' flag means that we are only intereseted in the snippet
% timestamps, not the actual snippet waveforms in this example.
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'snips', 'scalars'}, 'SORTNAME', SORTID, 'CHANNEL', CHANNEL, 'NODATA', 1);

%% Use TDTfilter to extract data around our epoc event
% Using the 'TIME' parameter extracts data only from the time range around
% our epoc event.
raster_data = TDTfilter(data{1,1}, REF_EPOC, 'TIME', TRANGE);

%%
% Adding the 'TIMEREF' flag makes all of the timestamps relative to the
% epoc event, which is ideal for generating histograms.
%hist_data = TDTfilter(data, 'Freq', 'VALUES', [100, 300]);filt by freq
hist_data = TDTfilter(data{1,1}, REF_EPOC, 'TIME', TRANGE, 'TIMEREF', 1);

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

%cdfidx=find((hist_data.epocs.Freq.data)>(cf-cf*0.2)&(hist_data.epocs.Freq.data)<(cf+cf*0.2));
hist_TS = hist_data.snips.(SNIP_STORE).ts;
subplot(2,1,1);
NBINS = floor(numel(hist_TS)/10);
hist(hist_TS, NBINS);
N = hist(hist_TS, NBINS); hold on;

axis tight; axis square;
set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
ylabel('Count','FontSize',16);
title({'Peristimulus', sprintf('Channel %d, n = %d trials', CHANNEL, num_trials)});

% % Draw a vertical line at t=0.
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
subplot(2,1,2);
plot(all_X, all_Y, '.', 'MarkerEdgeColor','k', 'MarkerSize',15); hold on;
axis tight; axis square;
set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
xlabel('Trial Window, s','FontSize',16);
ylabel('Trial Number','FontSize',16);
title({'Raster', sprintf('Channel %d, n = %d trials', CHANNEL, num_trials)});
% 
% % Draw a vertical line at t=0.
line([0 0], [0, trial+1], 'Color','r', 'LineStyle','-', 'LineWidth', 3);


%% Extracting spike times at a single lvl for CDF/latency calc.

 % get the data out and sort the responses by tone level
%[freq levl nlvl nLvls Responses TT figure1] = DataExtraction(MyTank,MyBlock);

% hopefully a better way to get levels, DataExtraction is old
levl = data.epocs.Levl.data;
[lvls,alpha]=TLevlFindNOrder(levl);

% % set the single lvl at which we want to calc. latency
%cdflvl=dBtoSPL(threshV,1); % convert it to dB
levl_cdf = levl;
threshV=max(levl);
% levl_cdf = sym(levl(1,:),'d'); % convert to sym because we want to round
% levl_cdf = round(levl_cdf*dec); % put in decimal for rounding at beginning of script
% levl_cdf=double(levl_cdf/dec); % convert back to double so it's not a useless imaginary number


% % find every trial where lvl == threshV and extract those trials
cdfSpikeTimes={};
tol = 0.00002; %tolerance for comparing the two numbers
for tr_ct=1:length(levl) % looping through the trials
    if ismembertol((levl(tr_ct,1)),threshV,tol)
       cdfSpikeTimes{tr_ct,1}=all_TS{tr_ct,1};
    end
end
% now cdfSpikeTimes has the spike times from the right trials

%removing empties
%cdfSpikeTimes(cellfun('isempty',cdfSpikeTimes))=[];

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
 spontR=(spont*4)/1000; % put it in terms of spikes per second
 
% [h,stats]=cdfplot(sortedtimes);
% xvals=get(h,'Xdata');
% yvals=get(h,'Ydata');
 
 
%% Cumulative spike fxn
% count spikes in this  time range rel. to stimulus
pre=-0.25;
post=0.2;

% making a csf for every trial
spkcum=[];
cdfbin = ((pre*1000):(post*1000))/1000;
for CDFcount = 1:length(cdfSpikeTimes)
    for BINcount = 1:length(cdfbin)
        spkcum(CDFcount,BINcount) = sum((cdfSpikeTimes{CDFcount}<cdfbin(BINcount)));
    end
%     spontR=spontR+spontR;
end





% avg across bins to get one csf
avgspkcum=[];
for avgct=1:length(spkcum)
    avgspkcum(avgct,1)=mean(spkcum(:,avgct));
end

% cumulative spont activity estimate
cumspont=zeros(length(avgspkcum),1);
for cumspontct=1:length(avgspkcum)
    cumspont(cumspontct,1)=cumspont(cumspontct,1)+spontR;
end

% subtracting spont from cumulative spike functions (Rowland et al. 2007)
for ff=1:length(avgspkcum)
    avgspkcum(ff,1)=avgspkcum(ff,1)-cumspont(ff,1);
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
freq=data.epocs.Freq.data;
mode(freq)