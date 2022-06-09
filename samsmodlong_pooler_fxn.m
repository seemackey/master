%% Modulation Transfer function analysis
% MTF analysis/plotting written by Sam Hauser (~2015), pooling and ROC analysis written by
% Chase Mackey (2021)
%
% 
% this samsmodlong_pooler_looper loops through this function
% to sample single-trials and do the johnson et al. (2012, j neurophys) 
% "across cell" method

% 
%

%%
function [pc] = samsmodlong_pooler_fxn(paths,data,VS)
% set variables

stdidx=24; % comparison stim for rate
bmf=28; %bmf for rate
stdvs=24; %comparison stim for VS 
bmfvs=27; %bmf for vs
flipr=0; %flip the sign if you're using downward sloping part of MTF
flipvs=0; %flip sign for vector strength

% duration
int = 1; %interval for stepping through binned spikes (2 for 500 ms, 4, for 250 ms, etc.)
dur = 0.0625; %bin width in seconds (0.5, 0.25, 0.1, 0.05)

% make this 1 if you want to do ROC analysis
discanalysis=1;

% reps for permutation loop
perms = 20;

% size of population
subpop = 0; % make this 1 if you want to randomly select < all neurons
popsize = 34; % how many neurons?



%% this lets you randomly select neurons from the population
if subpop == 1
    
    % select neurons here
    allsize = size(paths(:,1));
    randunits = randi([1 allsize(1,1)],popsize,1);
    data_sub = {};
    
    %loop to make the subpopulation = random selection of the real population
    for randsel = 1:1:length(randunits)
        data_sub(randsel,1) = data(randunits(randsel,1),1); 
    end
    
    data=data_sub;
    
end

%
%% looping through all paths for rate
% loop through every neuron and sum responses

check = 1; 
perm=1;
unit=1;
popspike=[];
popspike_all=[];

for perm = 1:1:perms                            % permutation loop
    for unit = 1:1:size(data)                  % neuron loop
        
        % bin the spikes 
        [stimspc,sortedspikes] = mfSort(data{unit,1},dur,int);

        % loop through the MTF function to get responses at each mf
        [levelspikes,onefreq] = MTF(stimspc);

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
            
    end
    popspike_all(:,perm) = sum(popspike,2); %pass the summed resps to the final array
end



%% Neuronal discrimination analysis
if discanalysis ==1
    %% ROC analysis 

    [pc,pFA,pHit]=mfROC_pop(popspike_all,stdidx,bmf,VS); % pass data to ROC function
    pc(:,2) = onefreq(1,stdidx:1:bmf)'; % get the mod freqs

end
end
