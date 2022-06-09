%
%%% INPUTS
% Monkey - Monkey name as a string
% Date - Date as a string in YYMMDD format
% Block - block number as a string
% Loc - Location of recording as string 'CN' or 'IC'
%%% OUTPUTS
% Data(1,:) = neurometric data
% Data(2,:) = tone levels in dBSPL
% N = spike counts by trial
% neurothresh = threshold of neuron based on neurometric curve


function [levl N] = MonkeyNandLevel(Monkey, Date, Num, Loc,SaveData)
Num=num2str(Num);

Date1=Date(1:6);


MyTank = ['Z:\' Monkey ' data\NeuroBehavior\' Loc '\ex' Date1 '\tank' Date];
%MyTank = ['H:\' Monkey ' data\Behavior\ex' Date1 '\tank' Date]; %using this to look do plot for ram really quickly

%    MyTank = ['C:\Documents and Settings\TDT1\My Documents\Experiments\NeuroBehavior\' x '\ex' Date '\tank' Date];
MyBlock = ['~OurData-' Num];



[levl N ] = MapDataExtraction(MyTank,MyBlock, Loc, Monkey,SaveData); % pulls out data


end



function [levl N] = MapDataExtraction(MyTank,MyBlock, Loc, Monkey,SaveData)

% set tank and location and select block
% MyTank = ['F:\Alpha data\Neurobehavior\IC\ex' Date '\tank' Date];
% MyBlock = ['~OurData-' Num];
MyTank;
if exist(MyTank) ~= 7 %#ok<*EXIST>
    warndlg('The tank you are looking for does not exist, please check inputs', 'Warning!');
    levl=0;
    N=0;
    return
end
BlockDir=[MyTank,'\',MyBlock(2:end)];
if exist(BlockDir)~=7
    warndlg('The Block you have selected does not exist within the tank you have selected','Failure!');
    levl=0;
    N=0;
    return
end


%Open tank/block for reading
figurea = figure('Position',[0 100000 1 1]);
TT = actxcontrol('TTank.X');
invoke(TT,'ConnectServer','Local','Me');
invoke(TT,'OpenTank',MyTank,'R');
invoke(TT,'SelectBlock',MyBlock);
invoke(TT,'ResetFilters');
Responses = invoke(TT, 'GetEpocsV', 'Corr', 0, 0, 10000); %puts the 'Corr' variable into an x,1 matrix with lots of zeros
channel = 1;
invoke(TT,'SetGlobalV','Channel',channel);
freq = invoke(TT, 'GetEpocsV', 'Freq', 0, 0, 10000);
global freqout
freqout=rot90(freq(1,:));
levl = invoke(TT, 'GetEpocsV', 'Levl', 0, 0, 10000);
nlvl = invoke(TT, 'GetEpocsV', 'nlvl', 0, 0, 10000);
nLvls= invoke(TT, 'GetEpocsV', 'NseL',0, 0, 30000);
%  nlvl = invoke(TT, 'GetEpocsV', 'Levl', 0, 0, 10000);
%  levl = invoke(TT, 'GetEpocsV', 'nlvl', 0, 0, 10000);

% length(Responses)
% length(freq)
% length(levl)
% length(nlvl)
% unique(Responses(1,:))
part1=length(freq(1,:));
part2=length(unique(freq(1,:)));
BFused=mode(freq(1,:));

part4=min(freq(1,:));
warn1str=num2str(part1);
warn1str=['Only ',warn1str, ' Trials Recorded, Terminating Analysis'];
if part1<70
    warndlg(warn1str, 'Warning!');
    close(figurea);
    levl=0;
    N=0;
    return
elseif part2>3 && part4>0
    warndlg('Block selected is most likely Response Map', 'Warning!');
    close(figurea);
    levl=0;
    N=0;
    return
end

%tlvls = TLevlFindNOrder(levl)
tlvls = unique(levl(1,:)); % contains tone levels used

CheckIfIsMod=find(freq(1,:)==12321);
if isempty(CheckIfIsMod)==0
    ModNoise=1;
    freq(1,CheckIfIsMod)=median(freq(1,:));
else
    ModNoise=0;
end

ROCFreq=median(freq(1,:));  %the frequency for the ROC analysis
ROCNoise=median(nlvl(1,:)); %ROC Noise Level
if ROCNoise>0
    v0=.79/(10^(82/20));
    ROCNoise=round(20*log10((ROCNoise)/v0));
elseif  ROCNoise==0
    ROCNoise=24;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psychometric function calc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ProbabilityCorrect=MonkeyTankPlotter2(freq,levl,nlvl,Responses,TT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x = 'Freq>=0';
TT.SetFilterWithDescEx(x);
trs = TT.GetValidTimeRangesV; %gets start and end time of each epoch 
TT.SetFilterTolerance(0.0001); %%If this code breaks, look at the returns for these functions to see if they are 1/0
TT.SetEpocTimeFilterV('Levl',0,.2);
TT.SetFilterWithDescEx(x);


maxtime = 1/325;
[~,y] = size(levl);
[~,y1] = size(trs);
if y > 400 % deletes extraneous trials
    y = 400;
    freq = freq(:,1:400);
    levl = levl(:,1:400);
    nlvl = nlvl(:,1:400);
end

if y>y1 % deletes extraneous trials
    y = y1;
    freq = freq(:,1:y1);
    levl = levl(:,1:y1);
    nlvl = nlvl(:,1:y1);
end

r = rem(y,10);
if r>0 % deletes extraneous trialsclc
    y = y-r;
    freq = freq(:,1:y);
    levl = levl(:,1:y);
    nlvl = nlvl(:,1:y);
end

N = zeros(1,y);
M=N;
DeleteCount = zeros(1,y);
SnipTimes = zeros(y,20);
InstRates = zeros(y,20);


%%%%%%%%%%%%%%%% CHANGE DURATION OF ANALYSIS WINDOW HERE %%%%%%%%%%%%%%
% dur=0.00325;
% start=0.05; % starting at latency, calculated in Raster_PSTH_latency
% 
% trs(1,:)=trs(1,:)+start; %changing start of trial 
% %trs(2,:)=(trs(1,:)+start+dur);
% trs(2,:)=trs(2,:)-(0.2-(dur+start)); %changing end of trial
% windowsize=trs(2,1)-trs(1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 
for l = 1:y %y is the trial count
    N(l) = TT.ReadEventsV(100000,'Snip',channel,0,trs(1,l),trs(2,l),'FILTERED'); 
    
    SnipTimes(l,1:N(l)) = TT.ParseEvInfoV (0,N(l),6);
    m=2;
    sniplength = N(l);
    while m <= sniplength % deletes spikes not neurologically possible
        datan = SnipTimes(l,m)-SnipTimes(l,m-1);
        if datan<maxtime
            SnipTimes(l,1:N(l))=[SnipTimes(l,1:m-1) SnipTimes(l,m+1:N(l)) 0];
            DeleteCount(l) = DeleteCount(l) + 1;
        else
            m = m+1;
        end
        sniplength = nonzerolength(SnipTimes(l,:));
    end
    InstRate = zeros(1,sniplength);
    InstRate(1) = 1/SnipTimes(l,1);
    for i = sniplength:-1:2
        InstRate(i) = 1/(SnipTimes(l,i)-SnipTimes(l,i-1));
    end
    InstRates(l,1:sniplength) = InstRate;
end


for l = 1:y %y is the trial count
    M(l) = TT.ReadEventsV(100000,'eNeu',channel,0,trs(1,l),trs(2,l),'FILTERED');
    
    SnipTimes(l,1:M(l)) = TT.ParseEvInfoV (0,M(l),6);
    
    m=2;
    sniplength = M(l);
    while m <= sniplength % deletes spikes not neurologically possible
        datan = SnipTimes(l,m)-SnipTimes(l,m-1);
        if datan<maxtime
            SnipTimes(l,1:M(l))=[SnipTimes(l,1:m-1) SnipTimes(l,m+1:M(l)) 0];
            DeleteCount(l) = DeleteCount(l) + 1;
        else
            m = m+1;
        end
        sniplength = nonzerolength(SnipTimes(l,:));
    end
    InstRate = zeros(1,sniplength);
    InstRate(1) = 1/SnipTimes(l,1);
    for i = sniplength:-1:2
        InstRate(i) = 1/(SnipTimes(l,i)-SnipTimes(l,i-1));
    end
    InstRates(l,1:sniplength) = InstRate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% LATENCY CALC. %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % set the single lvl at which we want to calc. latency 
% threshV=0.00353;
% cdflvl=dBtoSPL(threshV,1); % convert it to dB
% dec=1e5; % round to 5th decimal
% levl_cdf = sym(levl(1,:),'d'); % convert to sym because we want to round
% levl_cdf = round(levl_cdf * dec);
% levl_cdf=double(levl_cdf/dec); % convert back to double so it's not useless
% 
% 
% % % find every trial where lvl == threshV
% trs4cdf_idx=[];
% for tr_ct=1:length(levl) % y was previously defined as the num of trials
%     if levl_cdf(1,tr_ct)==threshV
%        trs4cdf_idx(tr_ct,1)=1;
%     end
% end
% % 
% % make a SnipTimes that has baseline activity in it (i.e. 250 ms before)
% start=-0.25; % start baseline calc. here
% % 
% trs_cdf=trs; %trs was previously defined as the trial start/end times
% trs_cdf(1,:)=trs_cdf(1,:)+start; %changing start of trial so we can look at baseline
% 
% 
% % get snip times with baseline rate incorporated 
% cdfSnipTimes=[];
% for snip_ct = 1:y %y is the trial count
%     tdt_extract(snip_ct) = TT.ReadEventsV(100000,'eNeu',channel,0,trs_cdf(1,snip_ct),trs_cdf(2,snip_ct),'FILTERED');
%     
%     cdfSnipTimes(snip_ct,1:tdt_extract(snip_ct)) = TT.ParseEvInfoV (0,tdt_extract(snip_ct),6);
%     
%     min_cdf=2;
%     sniplength_cdf = tdt_extract(snip_ct);
%     while min_cdf <= sniplength_cdf % deletes spikes not neurologically possible
%         datan_cdf = cdfSnipTimes(snip_ct,min_cdf)-cdfSnipTimes(snip_ct,min_cdf-1);
%         if datan_cdf<maxtime
%             cdfSnipTimes(snip_ct,1:tdt_extract(snip_ct))=[cdfSnipTimes(snip_ct,1:min_cdf-1) cdfSnipTimes(snip_ct,min_cdf+1:tdt_extract(snip_ct)) 0];
%             DeleteCount(snip_ct) = DeleteCount(snip_ct) + 1;
%         else
%             min_cdf = min_cdf+1;
%         end
%         sniplength_cdf = nonzerolength(cdfSnipTimes(snip_ct,:));
%     end
% end
% 
% cdfSnipTimes(cdfSnipTimes==0)=nan;
% % 
% % 
% % % go through the new SnipTimes and make a cdf of every response where lvl==threshV
% % cdf_yvals={};
% % cdf_xvals={};
% z=1;
% for cdfcount=1:length(trs4cdf_idx)
%     if trs4cdf_idx(cdfcount,1)==1
%      [h,stats]=cdfplot(cdfSnipTimes(trs4cdf_idx(cdfcount,:)));
%      cdf_yvals{cdfcount,z}=get(h,'YData');
%      cdf_xvals{cdfcount,z}=get(h,'XData');
%     end
%     z=z+1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% end of LATENCY CALC. %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

eNeuSnips=length(unique(M));
SnipSnips=length(unique(N));
if(eNeuSnips>SnipSnips)
    N=M; % puts all the spike counts from each trial into N, why??
end


% Chan1Snip = TT.ReadWavesOnTimeRangeV('Snip',1);
% N2=zeros(1,length(freq));
%close(figurea);
% length(freq(1,:))
% for i=1:length(freq(1,:))
%     N2(i)=length(unique(Chan1Snip(:,i)));
% end
% N=N2

Sorted = SortingByTime(levl, Responses);

m=length(Sorted(1,:));
for i=1:m-1
    Sorted(4,i)=N(i); % puts # spikes in each trial in the 4th row of "Sorted"
end


%%%%%%%%%%%%%%%%%%%
%                   Data has been collected at this point from the Data
%                   Tank
%
%%%%%%%%%%%%%%%%%%%
a=length(levl);
for i=1:a
    levl(2,i)=N(i);
end
b=length(tlvls);
ToneLevels=tlvls;
LengthTlvls=b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%converts voltage to dbSPL; probably dumb to put this here
% dblist = zeros(1,length(tlvls));
% conv = .79/(10^(82/20));
% 
% for i = 1:length(tlvls)
%     dblist(i) = 20*log10((tlvls(i))/conv);
% end
% dblist
% 
% approxthreshlvl=[7.89999976404943e-05];
% for trial_ct=1:length(y)
%     if Sorted(1,trial_ct) == approxthreshlvl
%         
%         

%%%%%%%%%%%%%%%%%%
%                   Build an array with Tone levels as the first column,
%                   assuming that the tone is repeated 30 times max. If
%                   this number is greater, make length the number of
%                   repeats +1
%%%%%%%%%%%%%%%%%%
for i=2:35
    tlvls(i,:)=NaN;
end


        
%%%%%%%%%%%%%%%%%%
%                   The following for loop sorts spike counts by tone from
%                   the sorted array into two arrays, tlvlsCorrect which
%                   corresponds to a correct pull from the monkey, and
%                   tlvlsIncorrect which contains all spike data for a no
%                   pull situation
%%%%%%%%%%%%%%%%%%
tlvlsCorrect=tlvls;
tlvlsIncorrect=tlvls;
for i=1:b
    count1=2;
    count2=2;
    for j=1:a
        if Sorted(1,j)==tlvls(1,i)
            if Sorted(3,j)==1
                tlvlsCorrect(count1,i)=Sorted(4,j);
                count1=count1+1;
            elseif Sorted(3,j)==2
                tlvlsIncorrect(count2,i)=Sorted(4,j);
                count2=count2+1;
            end
        end
    end
end


NANflag=0;
for i=1:length(tlvlsCorrect(1,:))
    Check1=tlvlsCorrect(2,i);
    Check2=tlvlsIncorrect(2,i);
    if isnan(Check1)==1 && isnan(Check2)==1
        NANflag=1;
        RowIndex=i;
    end
end
if NANflag==1
    tlvlsCorrect(:,RowIndex)=[];
    tlvlsIncorrect(:,RowIndex)=[];
end

AllDIs=zeros(1,100);
for zelda=1:100
    AllDIs(zelda)=MonkeyBootstrapDI(tlvlsCorrect,tlvlsIncorrect);
end

BootstrapDI=mean(AllDIs);
BootstrapDIstd=std(AllDIs);

[DI, WeightedDI]=FindDI(tlvlsCorrect,tlvlsIncorrect);

%%%%%%%%%%%%%%%%%%
%                   Sorting spikes into columns by tone level array at this
%                   point. Row one is tone level, each column will have the
%                   corresponding spike counts for each trial at that tone
%                   level
%
%%%%%%%%%%%%%%%%%%
for i=1:b
    c=2;
    for j=1:a
        if tlvls(1,i)==levl(1,j)
            tlvls(c,i)=levl(2,j);
            c=c+1;
        end
    end
end
if NANflag==1
    LengthTlvls=LengthTlvls-1;
    b=b-1;
    tlvls(:,RowIndex)=[];
end


%%%%%%%%%%%%%%%%%%
%                   Now that spikes are sorted, we want to make some
%                   parameters before we build an array for plotting a
%                   "histogram" these include the width histogram plots.
%                   For a smooth looking histogram, change the StepSizePlot
%                   value to a greater value. The default should be 5. The
%                   StepSize constant is the actual value used for ROC
%                   analysis. Currently it is set as low as possible to
%                   pool spikes  as accurately as possible.
%
%%%%%%%%%%%%%%%%%%
StepSize=1; %spike seperator
LowestSpikes=round((min(N)/10)-1)*10; %lowest spike count in the array N, which at this point has spikes counted
LowestPlotSpikes=round((min(N)/10)-1)*10;
StepSizePlot=5;%ceil((max(N)-min(N))/10);
HighestSpikes=round((max(N)/10)+.5)*10;
CountArray=(LowestSpikes:StepSize:HighestSpikes);
PlotArray=(LowestPlotSpikes:StepSizePlot:HighestSpikes);
LengthPlotArray=length(PlotArray);
LengthCountArray=length(CountArray);

%%%%%%%%%%%%%%%%%%
%                   Array for plotting histogram is constructed here with
%                   step sizes set from above
%%%%%%%%%%%%%%%%%%
for i=2:(b+1)
    PlotArray(i,:)=0;
end

for i=1:b
    for j=2:31
        for k=1:LengthPlotArray-1
            if tlvls(j,i)>=PlotArray(1,k) && tlvls(j,i)<PlotArray(1,k+1)
                PlotArray(i+1,k)=PlotArray(i+1,k)+1;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%
%                       Make array for ROC analysis with default step size
%                       of 1 as set above
%%%%%%%%%%%%%%%%%%%%
for i=2:(b+1)
    CountArray(i,:)=0;
end

for i=1:b
    for j=2:31
        for k=1:LengthCountArray-1
            if tlvls(j,i)>=CountArray(1,k) && tlvls(j,i)<CountArray(1,k+1)
                CountArray(i+1,k)=CountArray(i+1,k)+1;
            end
        end
    end
end
a=length(CountArray);



[RSmean, RSstd]= MonkeyBStrapResamp(CountArray,tlvls(1,:));
ROCarray=zeros(LengthTlvls,a);
%%%%%%%%%%%%%%%%%%%%
%                       At this point data is ordered in an array by bin numbers
%                       with # spikes in each bin into CountArray. The next
%                       for loop goes row by row, element by element
%                       finding the proportion of spikes to the left over
%                       total number of spikes in a row. Plotting any row
%                       against row 1 will produce an ROC.
%
%%%%%%%%%%%%%%%%%%%%
for j=1:LengthTlvls
    for i=1:a
        d=sum(CountArray(j+1,1:i));
        e=sum(CountArray(j+1,:));
        ROCarray(j,i)=1-(d/(e+eps));
    end
end

%%%%%%%%%%%%%%%%%%%%
%                       Finding area under the curve using rectangle and
%                       triangle method.
%%%%%%%%%%%%%%%%%%%%
lengthROC=length(ROCarray(1,:));
heightROC=length(ROCarray(:,1));
sums=zeros(heightROC-1,2);

for j=2:heightROC
    for i=(lengthROC):-1:2
        X1=ROCarray(1,i);
        X2=ROCarray(1,i-1);
        Y2=ROCarray(j,i-1);
        Y1=ROCarray(j,i);
        deltaX=X2-X1;
        deltaY=Y2-Y1;
        sums(j-1,2)=sums(j-1,2)+Y1*deltaX+deltaX*deltaY/2;
        sums(j-1,1)=ToneLevels(1,j);
    end
end

%Sums is area under the curve. Also corresponds to P(C).
sums=rot90(sums,3);
sums=fliplr(sums); %being lazy and orienting sums the same way as ProbabilityCorrect for neurometric data
alpha=length(sums(1,:)); % not sure why but having problem with list... chaning sums to tlvls


v0=.79/(10^(82/20));
for i = 1:alpha
    sums(1,i) = 20*log10((sums(1,i))/v0);
end

z=length(DI(1,:));
for i=1:z
    DI(1,i)=20*log10((tlvls(1,i))/v0);
end





sums=flipud(sums);
if SaveData==1
    
    ROCDataOutput4(ROCarray,MyTank,MyBlock,tlvls(1,:));
end

[BehaviorThresh,BehaviorSlope]=FindThreshold(ProbabilityCorrect(2,2:end),ProbabilityCorrect(1,2:end));
[BehavThreshFit,~,BehavXFit,BehavYFit]=FindFitThreshold(ProbabilityCorrect(2,2:end),ProbabilityCorrect(1,2:end));
[NeuroThresh,NeuroSlope]=FindThreshold(sums(2,:), sums(1,:));
[NeuroThreshFit,~,NeuroXFit,NeuroYFit]=FindFitThreshold(sums(2,:), sums(1,:));




[SpikeSlope,SpikeR, SpikedBSl, SpikedBR]=MonkeyFindSpikeSlopes(Sorted,NeuroThresh-abs(RSstd));

if SaveData==1
    
    ROCDataOutput1(PlotArray,MyTank,MyBlock,tlvls,Sorted,NeuroThresh-abs(RSstd));
    ROCDataOutput2(sums, MyTank,MyBlock);
    ROCDataOutput3(ProbabilityCorrect(:,2:end), MyTank,MyBlock);
    
end

RSBstd=0;

% if SaveData==1;
%     SaveThreshData(MyTank,MyBlock,NeuroThresh(1),BehaviorThresh(1),...
%         NeuroSlope,BehaviorSlope, WeightedDI, Loc, ROCNoise, Monkey,ModNoise,...
%         RSstd,SpikeSlope,SpikeR,SpikedBSl, SpikedBR, BFused,...
%         BehavThreshFit(1),NeuroThreshFit(1),RSBstd,BootstrapDI,BootstrapDIstd)
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz = get(0,'ScreenSize'); %ScreenSize is a four-element vector: [left, bottom, width, height]:
figure('Position',[50 scrsz(4)/4 scrsz(3)/1.5 (scrsz(4)/1.5-75)])
subplot(2,2,1),
hold on
if length(PlotArray(:,1))>=2
    plot(PlotArray(1,:),PlotArray(2,:),'k')
end
if length(PlotArray(:,1))>=5
    plot(PlotArray(1,:),PlotArray(5,:),'r')
end
if length(PlotArray(:,1))>=10
    plot(PlotArray(1,:),PlotArray(10,:),'c')
end
if length(PlotArray(:,1))>=14
    plot(PlotArray(1,:),PlotArray(14,:),'g')
end
ylabel('Occurrences');
xlabel('Spikes');
hold off


subplot(2,2,2)
hold on
if length(ROCarray(:,1))>=1
    plot(ROCarray(1,:),ROCarray(1,:),'k')
end
if length(ROCarray(:,1))>=4
    plot(ROCarray(1,:),ROCarray(4,:),'r')
end
if length(ROCarray(:,1,:))>=9
    plot(ROCarray(1,:),ROCarray(9,:),'c')
end
if length(ROCarray(:,1))>=13
    plot(ROCarray(1,:),ROCarray(13,:),'g')
end
ylabel('HR');
xlabel('FA');
hold off


subplot(2,2,3)
hold on
plot(sums(2,:),sums(1,:),'+k')
plot(NeuroThresh(1),NeuroThresh(2),'*k')
plot(NeuroXFit,NeuroYFit,'--k')
plot(NeuroThreshFit(1),NeuroThreshFit(2),'*b')

plot(ProbabilityCorrect(2,2:end),ProbabilityCorrect(1,2:end), '+r')
plot(BehaviorThresh(1),BehaviorThresh(2),'*r')
plot(BehavXFit,BehavYFit,'--r')
plot(BehavThreshFit(1),BehavThreshFit(2),'*g')
ROCNoise=num2str(ROCNoise);
ROCFreq =num2str(ROCFreq);
if ModNoise==0;
    a=['P(C) at ' ROCFreq 'Hz in ' ROCNoise 'dB noise'];
elseif ModNoise==1;
    a=['P(C) at ' ROCFreq 'Hz in ' ROCNoise 'dB Modulated noise'];
end
title(a ,'FontWeight','bold')
ylabel('Percent');
xlabel('Tone level (dB), (neuro in black)');
legend off
hold off

subplot(2,2,4)
DI(1,1)=-99;
hold on
plot(DI(1,:),DI(2,:), '*k')
plot(0,WeightedDI,'*r');
hold off
WeightedDIStr=num2str(WeightedDI);
a='Discriminate Indices';
title(a ,...
    'FontWeight','bold')
ylabel('DI');
xlabel('Tone Level in dB');

text(10,WeightedDI,WeightedDIStr)

end

function [LevelAndResponses] = SortingByTime(levl, Responses)




a=length(levl);
b=length(Responses);
j=1;

for i=1:a
    while(levl(2,i)>=Responses(2,j))
        j=j+1;
        if j>b
            LevelAndResponses=levl;
            return
        end
        
    end
    levl(3,i)=Responses(1,j);
    levl(5,i)=Responses(2,j)-levl(2,i);
end
LevelAndResponses=levl;
end


function [DI,WeightedDI]=FindDI(Correct,Incorrect)


a=length(Correct(1,:));
b=length(Correct(:,1));
DI=zeros(3,a);
DI(1,:)=Correct(1,:);
Greater=zeros(1,a);
Equal=zeros(1,a);
Dimension=zeros(2,a);
% Correct
% Incorrect

Correct(1,:)=-Correct(1,:);        %Making the first row of correct and incorrect
Incorrect(1,:)=-Correct(1,:);      %inverted.



for i=1:a %a is height correct/incorrect
    for j=2:b %b is length correct/incorrect
        if isnan(Correct(j,i))==0
            Dimension(1,i)=Dimension(1,i)+1;
        end
        if  isnan(Incorrect(j,i))==0
            Dimension(2,i)=Dimension(2,i)+1;
        end
    end
end

for i=1:a
    j=2;
    for l=2:b
        Comp1=Correct(j,i);
        for k=2:b
            if Comp1 > Incorrect(k,i)
                Greater(1,i)=1+Greater(1,i);
            end
            if Comp1 == Incorrect(k,i)
                Equal(1,i)=Equal(1,i)+1;
            end
        end
        j=j+1;
    end
end

SumWeight=0;
for i=1:a
    DI(2,i)=(Greater(i)+.5*Equal(i))/(Dimension(1,i)*Dimension(2,i));
    if isnan(DI(2,i))==0;
        if i==1
            DI(3,i)=DI(2,i)*Dimension(1,i);
            SumWeight=Dimension(1,i)+SumWeight;
        else
            DI(3,i)=DI(2,i)*Dimension(2,i);
            SumWeight=Dimension(2,i)+SumWeight;
        end
    end
end

WeightedDI=sum(DI(3,:)/SumWeight);
end

function [ThreshInDB, slope]=FindThreshold(XDataIn, YDataIn)
threshold=.76;
a = length(YDataIn);
YDifference=zeros(1,a-1);


for i=2:a
    if YDataIn(i-1)<.76 && YDataIn(i)>=.76
        YDifference(i-1)=YDataIn(i)-YDataIn(i-1);
    else
        YDifference(i-1)=NaN;
    end
end


[Value,Indices]=max(YDifference);
slope=Value/(XDataIn(Indices+1)-XDataIn(Indices));
b=YDataIn(Indices)-slope*XDataIn(Indices);
xvalues=(XDataIn(Indices):.1:XDataIn(Indices+1));
yvalues=slope*xvalues+b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           To find the closest value to threshold, what you can do is subtract .76
%                           from all numbers, find the absolute value of all numbers, and then find
%                           the minimum index value in the array. This will find the closest
%                           intercept of a non continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,index]=min(abs(yvalues-threshold));
ThreshInDB(1)=(xvalues(index)*slope)/slope;
ThreshInDB(2)=(yvalues(index)*slope)/slope;



end
function [ThreshInDB,slope,x,y]=FindFitThreshold(XDataIn, YDataIn)
threshold=.76;
temp = min(XDataIn);
xtralvls=rot90(XDataIn)
XDataIn=rot90(XDataIn)-temp;
YDataIn=rot90(YDataIn)

try
    coeff=CreateWeibulFit(XDataIn,YDataIn);
    XDataIn=XDataIn+temp;
    a=coeff(1);
    b=coeff(2);
    %     c=coeff(3);
    d=coeff(3);
    x=(min(XDataIn):0.1:max(XDataIn)); % original step size ((max(XDataIn)-min(XDataIn))/500)
    y=WeibulCurveFitFunction(x,a,b,d,temp);  %createFit3 is thus far the best fit over createFit.
    YDifference=zeros(1,length(y)-1);
    for i=2:length(y)
        if y(i-1)<=threshold && y(i)>threshold
            YDifference(i-1)=y(i)-y(i-1);
        else
            YDifference(i-1)=NaN;
        end
    end
    
    [Value,Indices]=max(YDifference);
    slope=Value/(x(Indices+1)-x(Indices));
    [~,index]=min(abs(y-threshold));
    ThreshInDB(1)=x(index);
    ThreshInDB(2)=y(index);
catch
    ThreshInDB=[NaN,NaN];
    slope=NaN;
    x=NaN;
    y=NaN;
end
end
