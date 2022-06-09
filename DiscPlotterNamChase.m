%function [] = DiscPlotter(Monkey, Date, Num)
% this function takes a block of discrimination data and does our typical analysis of
% accuracy and speed 
% Chase Mackey, Sept. 2020
%edited by Namrata 09/14/2020
tic
clear all;close all;
%
Monkey='Dario';
Date='210611d';
Num='3';
Num=num2str(Num);
Date1=Date(1:6);
% want to import the stimuli? it takes a while. 1 for yes.
gimmestimmies=0;
%/Volumes/ramlab2018/Isildur data/Behavior/ex200818
MyTank = ['/Volumes/ramlab2018/' Monkey ' data\Behavior\ex' Date1 '\tank' Date];
MyBlock = ['\OurData-' Num];
path = ['/Volumes/ramlab2018/' Monkey ' data/Behavior/ex' Date1 '/tank' Date '/OurData-' Num];

%% Data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = TDTbin2mat(path);

% put everything into a matrix
% but first check if trial num = response num
Trs = length(Data.epocs.Freq.data); %trials
Resps = length(Data.epocs.Corr.data); %Responses
if isequal (Trs,Resps)
    thingswereequal=1;
    All(1,:) = Data.epocs.Freq.data'; %trials
    All(2,:) = Data.epocs.CaTr.data'; % 0 is a catch trial, 1 is a signal trial
    All(3,:) = Data.epocs.S__2.data'; % key parameter of second stimulus (e.g. mod freq)
    All(4,:) = Data.epocs.Corr.data'; %Responses
    All(5,:) = Data.epocs.Corr.onset'; % Response time
    All(6,:) = Data.epocs.S__2.onset'; %S2 onset time
    %All(7,:) = Data.epocs.Erly.data'; % early releases, this might not be working in the circuit yet.
else
    % delete last trial to compensate (triple checked that this is right)
    Data.epocs.Corr.data(1,:)=[];
    Data.epocs.Corr.onset(1,:)=[]; % Response time
    Data.epocs.Corr.offset(1,:)=[];
    %Data.epocs.Erly.data(1,:)=[];
    %Data.epocs.Erly.onset(1,:)=[]; % Response time
    %Data.epocs.Erly.offset(1,:)=[];
    % then get the data in one matrix now that things are equal
    All(1,:) = Data.epocs.Freq.data'; %trials
    All(2,:) = Data.epocs.CaTr.data'; % 0 is a catch trial, 1 is a signal trial
    All(3,:) = Data.epocs.S__2.data'; % key parameter of second stimulus (e.g. mod freq)
    All(4,:) = Data.epocs.Corr.data'; %Responses
    All(5,:) = Data.epocs.Corr.onset'; % Response time
    All(6,:) = Data.epocs.S__2.onset'; %S2 onset time
    %All(7,:) = Data.epocs.Erly.data'; % early releases, this might not be working in the circuit yet.
end
if gimmestimmies==1
    Stim(1,:) = Data.streams.stim.data'; %stimuli
    time=1:1:length(Stim)'; %time variable in ms for stimuli
end
    
    % find the duration of the stimulus
    Durs=Data.epocs.Freq.offset-Data.epocs.Freq.onset;
    Duration=median(Durs); % only use this as your dur if the dur is constant
    
%%
% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% false alarm rate calc.
ct = length(find(All(2,:)==0));   %%EASIER WAY
 
%  counts up # of false alarms
fa=0;
for FA_count=1:length(All)
    if All(4,FA_count)==1 && All(2,FA_count)==0
        fa=fa+1;
    else 
        continue
    end
end

% divide them both to get false alarm rate
FArate=fa/ct

%% calculate hit rate

Responses = All(4,:);
% find the stim diff for each trial
for i = 1:length(All(1,:))
diff(1,i) = (All(3,i) - All(1,i));
end

uniqueDif1 = unique(diff); %list of stimuli
uniqueDiff = uniqueDif1(1,2:end);   %%REMOVE catch trials
rr=[];
% counts up hits
for j = 1:length(uniqueDiff(1,:))
    ToneDiff(1,j) = uniqueDiff(1,j);
    stimtrialidx = find(diff ==(uniqueDiff(1,j))); %get index of every stim trial (not catch trials)
    l = 1;
    for k = 1:length(stimtrialidx(1,:)) %loop through stim trials
        rr(1,l) = Responses(1,stimtrialidx(1,k)); % store in temp variable rr
        l = l+1;
        
    end
    hitidx = find(rr == 1); % index of each stim trial where they released
    TF = isempty(hitidx);
    if TF == 1
         ProbCor(1,j)=0
    else
    ProbCor(1,j) = length(hitidx(1,:)) / length(stimtrialidx(1,:)); %num hits/num occurrence
    end
    rr = [];  %%%RESET RR
end 
  


%% calc dprime

%change the values because dprime computes infinity for 0 and 1
for reduce_ct = 1:length(ProbCor)
    if ProbCor(1,reduce_ct) == 1 
        ProbCor(1,reduce_ct) = 0.99;
    end
    if ProbCor(1,reduce_ct) == 0
        ProbCor(1,reduce_ct) = 0.01;
    end
end
dprime=norminv(ProbCor)-norminv(FArate);

% probability correct in 2AFC
ProbCor_2AFC=normcdf(dprime./2); % convert to prob correct


%% Reaction Time analysis


% RT calc. (lever release minus stim onset)
for RT_count = 1:length(All)
    if All(4,RT_count) == 1  %if resp was yes, doesn't care if it's a catch trial
    all_RTs(RT_count,1) = All(5,RT_count)-All(6,RT_count); % the RT
    all_RTs(RT_count,2) = All(3,RT_count)-All(1,RT_count); % the stim diff
    end
end


%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ftsize=18;
% make the FA rate a line on the hit rate plot
FA_x_axis = ToneDiff(1,1):1:ToneDiff(end);
figure

% plots hit rate
subplot(2,2,1)
plot(ToneDiff,ProbCor,'o')
H=gca;
H.LineWidth=1.5; 
hold on
plot(FA_x_axis,FArate,'x')
ylabel('Hit Rate')
xlabel('Stim Diff (Hz)')
axis square

% plot d prime even though it's stupid
subplot(2,2,2)
plot(ToneDiff,ProbCor_2AFC,'o')
H=gca;
H.LineWidth=1.5; 
ylabel('P(C)')
xlabel('Stim Diff (Hz)')
axis square

% plot all RTs (incl. catch trials)
subplot(2,2,3)
scatter(all_RTs(:,2),all_RTs(:,1))
H=gca;
H.LineWidth=1.5; 
ylabel('Reaction Time (s)')
xlabel('Stim Diff (Hz)')
axis square




%% RT cdf plotting 

% sort the RTs by mod freq. so we can make a few cdfs
RTminidx=find(all_RTs(:,2)==8); % kinda think this should be hard coded
RTmididx=find(all_RTs(:,2)==16);
RTmaxidx=find(all_RTs(:,2)==64);

% plot cdfs of those RTs which fall near, slightly above, and way above
% thresh
subplot(2,2,4)
cdfplot(all_RTs(RTminidx,1))
H=gca;
H.LineWidth=1.5; 
hold on
cdfplot(all_RTs(RTmididx,1))
H=gca;
H.LineWidth=1.5; 
hold on
cdfplot(all_RTs(RTmaxidx,1))
H=gca;
H.LineWidth=1.5; 
ylabel('Cum. Prop.')
xlabel('RT (s)')
axis square
legend('8 Hz','16 Hz','64 Hz')
%saveas(figure(1),[pwd '\DiscData\' Monkey Date Num '.fig']);

% plot accuracy w fit
thresh = ToneDiff(1,1);
% createFitWB_disc(ToneDiff,dprime,thresh) % function that does the weibull fit
createFitWB_HRdisc(ToneDiff(1,1:end),ProbCor_2AFC(1,1:end),thresh)
createFitWB_HRdisc(ToneDiff(1,1:end),ProbCor(1,1:end),thresh)

Duration
toc
% RT VS LEVEL FIT THAT DOESNT WORK%%%%%%%%%%%%%%%%%%%%%%%%%
% x=all_RTs(:,2);
% y=all_RTs(:,1);
% p=polyfit(all_RTs(:,2),all_RTs(:,1),1);
% r = p(1).*(all_RTs(:,2)) + p(2);
% figure
% plot(all_RTs(:,2),r)
% RTslope=p(1)
% RTint=p(2)


%end