%% Modulation Transfer function analysis
% MTF analysis/plotting written by Sam Hauser (~2015), pooling and ROC analysis written by
% Chase Mackey (2021)
%
% 
% this is the version that loops through multiple blocks. trying to make a
% CN--IC pooling model - CM
% now i'm adding onto the loop version so it adds responses from each
% neuron instead of just pooling binned spikes - CM
%
clear all; 
close all;

%% set variables

%respNonPref = onefreq(1,19); % setting the comparison stimulus
stdidx=15; % comparison stim for rate
bmf=18; %bmf for rate
stdvs=24; %comparison stim for VS 
bmfvs=27; %bmf for vs

flipr=0; %flip the sign if you're using downward sloping part of MTF
flipvs=0; %flip sign for vector strength

% duration
int = 2; %interval for stepping through binned spikes (2 for 500 ms, 4, for 250 ms, etc.)
dur = 1/int; %bin width in seconds (0.5, 0.25, 0.1, 0.05)

% make this 1 if you want to do ROC analysis
discanalysis=1;

% reps for permutation loop
perms = 20;

% size of population
subpop = 1; % make this 1 if you want to randomly select < all neurons
popsize = 32; % how many neurons?

%% file paths

paths = {
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181207/tank181207c/OurData-4';	% high mf units
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex180523/tank180523c/OurData-5';	
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181016/tank181016c/OurData-23';	
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151218/tank151218c/OurData-13';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150611/tank150611c/OurData-7';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150612/tank150612c/OurData-15';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-2';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-3';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-11';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-17';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150724/tank150724c/OurData-8';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150807/tank150807c/OurData-5';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151125/tank151125c/OurData-7';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151125/tank151125c/OurData-11';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151218/tank151218c/OurData-13';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151106/tank151106c/OurData-18';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex180523/tank180523c/OurData-5';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181005/tank181005c/OurData-22';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181009/tank181009c/OurData-13';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181012/tank181012c/OurData-16';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181016/tank181016c/OurData-23';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181018/tank181018c/OurData-31';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181113/tank181113c/OurData-14';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-9';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-18';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-22';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181129/tank181129c/OurData-5';
        '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181130/tank181130c/OurData-6';
         '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-13';
         '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-17';
         '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-5';
         '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-10';
         '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-19';
         '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-23';
         '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181207/tank181207c/OurData-4';
        };
    
% import
for imp = 1:1:size(paths)
    data{imp,1} = TDTbin2mat(paths{imp,1});
end

%% this lets you randomly select neurons from the population
if subpop == 1
    % select neurons here
    allsize = size(paths(:,1));
    randunits = randi([1 allsize(1,1)],popsize,1);
    data_sub = {};
    
    %loop to make the subpopulation = random selection of the real population
    for randsel = 1:1:length(randunits)
        data_sub(randsel,1) = data(randunits(randsel,1),1); %not indexing this correctly yet
    end
    
    data=data_sub;
    
end
%


%% looping through all paths
% loop through every neuron and sum responses
check = 1; 
for perm = 1:1:perms                            % permutation loop
    for unit = 1:1:size(data)                  % neuron loop

        % bin the spikes so we can make an MTF
        [binnedspc,stimspc,sortedspikes] = mfSort(data{unit,1},0.5,2);

        % loop through the MTF function to get responses at each mf
        [levelspikes,levelspikesVS,onefreq] = MTF(binnedspc,0.5,2,stimspc,sortedspikes);

        if check == 1
            for mf = 1:1:length(levelspikes)
            r = randi([1 4],1,1); % select a random trial, we had 4 reps
            popspike(mf,unit) = levelspikes(mf,r);
            end
            check = check+1;
        else 
            for mf = 1:1:length(levelspikes)
            r = randi([1 4],1,1); 
            popspike(mf,unit) = levelspikes(mf,r); %concatenate with existing array
            end
        end
            
        % add resps from each single unit's MTF 
        % resps are chosen randomly from each trial and from each neuron, resulting in
        % a distribution of 20 summed resps.

        % ROC on population MTF at a few MFs. 
    end
    popspike_all(:,perm) = sum(popspike,2); %pass the summed resps to the final array
end


%% Neuronal discrimination analysis
if discanalysis ==1
    %% ROC analysis on spike count

    [pc,FA,HR]=mfROC_pop(popspike_all,stdidx,bmf,0); % pass data to ROC function
    pc(:,2) = onefreq(1,stdidx:bmf)'; % get the mod freqs


    %% ROC analysis on vector strength 

%      [PCvs,FAvs,HRvs]=mfROC(levelspikesVS(1:end-1,:),stdvs,bmfvs,1); % pass data to ROC function
%      PCvs(:,2) = onefreq(1,stdvs:bmfvs)'; % get the mod freqs

    %% plot neurometric functions with MTF data
    
    % fit with weibull to extract thresh and slope
    [ratethresh,ratethreshdelta,sloper]=createFitMF(pc(:,2),pc(:,1),pc(1,2),0,flipr);
%    
%    [VSthresh,VSthreshdelta,slopevs]=createFitMF(PCvs(:,2),PCvs(:,1),PCvs(1,2),0,flipvs); 
   
%    % printing these out
%    slopes=[sloper,slopevs];
%    DiscThresholds = [ratethresh,VSthresh,ratethreshdelta,VSthreshdelta, sloper,slopevs]
   
   % PLOTTING EVERYTHING TOGETHER
    figure
    subplot(1,2,1)
    semilogx(onefreq,popspike_all,'color',[0,0,0]+0.8);
    hold on
    semilogx(onefreq,mean(popspike_all,2),'b');
    title('Spike Count');
    xlabel('Modulation Frequency'); 
    ylabel('Firing Rate (Hz)')
%     subplot(2,2,2)
%     semilogx(onefreq, done);
%     title('Vector Strength'); 
%     xlabel('Modulation Frequency');
%     ylabel('Vector Strength'); 
    subplot(1,2,2)
    plot(pc(:,2),pc(:,1),'-o') % plot probability correct in a 2AFC
    title('Rate-based neurometric function'); 
    xlabel('Modulation Frequency');
    ylabel('Neurometric P(C)');
%     subplot(2,2,4)
%     plot(PCvs(:,2),PCvs(:,1),'-o') % plot probability correct in a 2AFC
%     title('VS-based neurometric function');
%     xlabel('Modulation Frequency');
%     ylabel('Neurometric P(C)');
%     
%     % saving these for manual curve fitting
%     x=[PCvs(:,2)]; % vector strength
%     y=[PCvs(:,1)];
%     
%     xr=[pc(:,2)]; % rate
%     yr=[pc(:,1)];
end

