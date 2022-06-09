% Response Block Plotter
% Written by Andy Hrnicek under guidance of Ram Ramachandran, PhD.
% March 2011
% 
% This function plots a rate level function for the neuron recorded in
% neurophysiology. It returns the neurometric data that can be plotted on
% it's own or with the psychometric function; as is done in
% 'plottogether.m'
%
% This code only runs the crudely sorted data done by the Sorter in the
% OpenWorkbench. It needs to be adapted to use the sorted data created in
% OpenSorter.
%
% If you want to do a rate level function, return the variable called
% spikemed and you can plot the median spike count by tone level versus
% dbSPL with the psychometric function in the 'plottogether.m' script
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


function [Data N neurothresh] = MonkeyRespBlockPlotter(Monkey, Date, Num, Loc,Type)
Num=num2str(Num);

Date1=Date(1:6);
%does not work for vowels because of line 53 (x = ['Freq>0'];). this is
%setting the filter for accessing TDT storage. I don't know of a good
%solution for fixing this to look at vowel data. It gets stuck because
%there is a data overload because the filter isn't specific enough if you
%try to use Freq and some other logical condition. The other variables are
%too broad to really use here though, see TDT documentation on how to use
%their functions for more specifics.

if Type==1
    MyTank = ['G:\' Monkey ' data\Behavior\ex' Date1 '\tank' Date];
else
    MyTank = ['G:\' Monkey ' data\Neurobehavior\ex' Date1 '\tank' Date];
end
MyBlock = ['~OurData-' Num];


[Data N neurothresh] = MapDataExtraction(MyTank,MyBlock); % pulls out data


end



function [Data N neurothresh] = MapDataExtraction(MyTank,MyBlock)

% set tank and location and select block
% MyTank = ['F:\Alpha data\Neurobehavior\IC\ex' Date '\tank' Date];
% MyBlock = ['~OurData-' Num];

%Open tank/block for reading
figurea = figure;
TT = actxcontrol('TTank.X');
invoke(TT,'ConnectServer','Local','Me');
invoke(TT,'OpenTank',MyTank,'R');
invoke(TT,'SelectBlock',MyBlock);
invoke(TT,'ResetFilters');

Responses = invoke(TT, 'GetEpocsV', 'Corr', 0, 0, 10000); %puts the 'Corr' variable into an x,1 matrix with lots of zeros

% Specify channel
% code = invoke(TT, 'GetEventCodes', 0);
% start = invoke(TT,'CurBlockStartTime')
% formstart = invoke(TT,'FancyTime',start ,'D/O/Y H:M:S.U W')
% notes = invoke(TT,'CurBlockNotes')
channel = 1;
invoke(TT,'SetGlobalV','Channel',channel);
freq = invoke(TT, 'GetEpocsV', 'Freq', 0, 0, 10000);
levl = invoke(TT, 'GetEpocsV', 'Levl', 0, 0, 10000);
nlvl = invoke(TT, 'GetEpocsV', 'nlvl', 0, 0, 10000);

tlvls = TLevlFindNOrder(levl);

x = 'Freq>0';
filt0 = TT.SetFilterWithDescEx(x);
trs = TT.GetValidTimeRangesV;

%% pulls out all snip counts to all trials and stores in 'N'
TT.SetFilterTolerance(0.0001); %%If this code breaks, look at the returns for these functions to see if they are 1/0
TT.SetEpocTimeFilterV('Levl',0,.2);
% sortname = [Date '-' Num ' sorted']; 
% TT.SetUseSortName(sortname);
TT.SetFilterWithDescEx(x);


maxtime = 1/325;
[x,y] = size(levl);
[x1,y1] = size(trs);
if y > 260 % deletes extraneous trials
    y = 260;
    freq = freq(:,1:260);
    levl = levl(:,1:260);
    nlvl = nlvl(:,1:260);
end
if y>y1 % deletes extraneous trials
    y = y1;
    freq = freq(:,1:y1);
    levl = levl(:,1:y1);
    nlvl = nlvl(:,1:y1);
end
r = rem(y,20);
if r>0 % deletes extraneous trials
    y = y-r;
    freq = freq(:,1:y);
    levl = levl(:,1:y);
    nlvl = nlvl(:,1:y);
end

N = zeros(1,y);
DeleteCount = zeros(1,y);
SnipTimes = zeros(y,20);
InstRates = zeros(y,20);
for l = 1:y
    N(l) = TT.ReadEventsV(100000,'Snip',channel,0,trs(1,l),trs(2,l),'FILTERED');
    SnipTimes(l,1:N(l)) = TT.ParseEvInfoV(0,N(l),6);
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
close(figurea);

% bin tone levels
% median N by level of tone
[ToneLevel alpha] = TLevlFindNOrder(levl);
spikemean = zeros(1,alpha); spikemed = zeros(1,alpha); dprime = zeros(1,alpha); spikestd = zeros(1,alpha);
for i = 1:alpha
    indeces = find(levl(1,:) == ToneLevel(i));
    spikes = [];
    for j = 1:length(indeces)
        spikes(j) = N(indeces(j));
    end
    
    % calculate median spikes, mean, and standard deviation for each tone
    % level as it goes through and stores them into their respective
    % variables
    spikemed(i) = median(spikes);
    spikemean(i) = mean(spikes);
    spikestd(i) = std(spikes);
    if i > 1
        % calculates the z-score for the values of spike means
        dprime(i) = ((spikemean(i) - spikemean(1))/(((spikestd(i)^2)+(spikestd(1)^2))^.5));  
    end
end

dbSPL = dBtoSPL(ToneLevel, alpha);


% if you want to do a rate level function, return the variable called
% spikemed and you can plot the median spike count by tone level versus
% dbSPL with the psychometric function in the 'plottogether.m' script
maybe = normcdf(dprime);
likely = find(maybe>.76); likely = likely(1);
m = (maybe(likely)-maybe(likely-1))/(dbSPL(likely)-dbSPL(likely-1));
b = maybe(likely)-m*(dbSPL(likely));
Data = [maybe; dbSPL];
neurothresh = (.76-b)/m;


fig = figure;
[AX,H1,H2] = plotyy(dbSPL,spikemed,dbSPL,maybe);
% sets axes and labesl them
set(get(AX(1),'Ylabel'),'String','RateLevel') 
set(get(AX(2),'Ylabel'),'String','Neurometric') 
xlabel('Tone Level (dBSPL)') 
title('Neurobehavior Plot') 
% creates line styles for neurometric and psychometric - blue is
% pscyhometric and black is neurometric
set(H1,'LineStyle','-','Color','b')
set(H2,'LineStyle','-','Color','k','Marker','*')
end

% finds the tone levels used in block and sorts them in ascending order
function [ToneLevel alpha] = TLevlFindNOrder(levl)

tonelevels=[]; %creates array
ToneLevel=sort(levl(1,:)); %this becomes an infinite loop without the sort...
while ~isempty(ToneLevel) %iterates through until ToneLevel is and empty array
tonelevels=[min(ToneLevel),tonelevels];
    while(ToneLevel(1)==tonelevels(1)) %deletes all instances of a value added to tonelevels
    ToneLevel(1)=[];
        if isempty(ToneLevel) %prevents error due to exceeding matrix dimensions
            break
        end
    end
end

for i=1:length(tonelevels) %reorders them for ease of use
    ToneLevel=[tonelevels(i) ToneLevel];
end

alpha = length(tonelevels);

end

% converts voltage levels to dBSPL for tone levels to be plotted 
function dbSPL = dBtoSPL(ToneLevels, alpha)

dbSPL = zeros(1,alpha);
v0 = .79/(10^(82/20));
for i = 1:alpha
    dbSPL(i) = 20*log10((ToneLevels(i))/v0);
end

end
