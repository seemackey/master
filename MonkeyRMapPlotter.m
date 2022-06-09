% RMapPlotter
% Written by Andy Hrnicek under guidance of P.I. Ram Ramachandran, PhD.
% Spring 2011
% 
% This code plots one level of a response map, or a response map one block
% at a time. It returns a 2xN array of frequencies and spike count of the
% neuron at that frequency. 
% 
% This code does not use the sorted data and only is set up to work on the
% crude sort done with the Sorter in the Open Workbench. It needs to be
% adapted to use re-sorted data. Look at 'DataExtraction3.m' for ideas on
% how to change it to look at filtered data. 
% 
% The red 'o' that is plotted is
% the average of spike counts for the catch trials and is a measure of
% spontaneous activity or baseline firing rate.
%
%%% INPUTS
% Monkey - Monkey name as a string
% Date - Date as a string in YYMMDD format
% Block - block number as a string
% x or Loc - Location of recording as string 'CN' or 'IC'
%%% OUTPUTS
% Data(1,:) = Frequencies, Data(2,:) = Spike Counts


function [OutData1, OutData2, LevelOut, flag, baseline] = MonkeyRMapPlotter(Monkey, Date, Num, Loc)

Num=num2str(Num);
Date1=Date(1:6);
if nargin < 4
%     MyTank = ['F:\' Monkey ' data backup\' Monkey ' data\Behavior\ex' Date '\tank' Date];
    MyTank = ['C:\My Documents\Experiments\Behavior\' Monkey '\ex' Date1 '\tank' Date];
    MyBlock = ['~OurData-' Num];
else
    MyTank = ['G:\' Monkey ' data\NeuroBehavior\' Loc '\ex' Date1 '\tank' Date];
%    MyTank = ['C:\Documents and Settings\TDT1\My Documents\Experiments\NeuroBehavior\' x '\ex' Date '\tank' Date];
    MyBlock = ['~OurData-' Num];
end

flag=0;
if exist(MyTank) ~= 7
    f=warndlg('The tank you are looking for does not exist, please check inputs', 'Warning!');
    waitfor(f);
    OutData1=NaN;
    OutData2=NaN;
    LevelOut=NaN;
    flag=1;
    return
end
BlockDir=[MyTank,'\',MyBlock(2:end)];
if exist(BlockDir)~=7
    f=warndlg('The Block you have selected does not exist within the tank you have selected','Failure!');
    waitfor(f);
    OutData1=NaN;
    OutData2=NaN;    
    LevelOut=NaN;
    flag=1;
    return
end


[Data OutData1 OutData2 LevelOut flag baseline] = MapDataExtraction(MyTank,MyBlock,flag); % pulls out data





end



function [Data OutData1 OutData2 LevelOut flag baseline] = MapDataExtraction(MyTank,MyBlock, flag)

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

% Specify channel
% code = invoke(TT, 'GetEventCodes', 0);
% start = invoke(TT,'CurBlockStartTime')
% formstart = invoke(TT,'FancyTime',start ,'D/O/Y H:M:S.U W')
% notes = invoke(TT,'CurBlockNotes')
channel = 1;
invoke(TT,'SetGlobalV','Channel',channel);
freq = invoke(TT, 'GetEpocsV', 'Freq', 0, 0, 10000);
freq(1,:);
levl = invoke(TT, 'GetEpocsV', 'Levl', 0, 0, 10000);
levl(1,:);

part1=length(unique(freq(1,:)));
if part1<5
    f=warndlg('Response map short or wrong block selected', 'Warning!');    
    close(figurea);
    waitfor(f);
    baseline=NaN;
    Data    =NaN;
    OutData1=NaN;
    OutData2=NaN;
    LevelOut=NaN;
    flag    =1;
    return
else flag=0;
end

LevelOut=mode(levl(1,:));
x = 'Freq>0';
filt0 = TT.SetFilterWithDescEx(x);
trs = TT.GetValidTimeRangesV;

%% pulls out all snip counts to all trials and stores in 'N'
TT.SetFilterTolerance(0.0001); %%If this code breaks, look at the returns for these functions to see if they are 1/0
TT.SetEpocTimeFilterV('Levl',0,.2);
% sortname = [Date '-' Num ' sorted']; 
% TT.SetUseSortName(sortname);
TT.SetFilterWithDescEx(x);
[x,y] = size(levl);
if y > 110
    y = 110;
end
N = zeros(1,y);
M= zeros(1,y);
for l = 1:y
    N(l) = TT.ReadEventsV(100000,'Snip',channel,0,trs(1,l),trs(2,l),'FILTERED');
end
for l = 1:y
    M(l) = TT.ReadEventsV(100000,'eNeu',channel,0,trs(1,l),trs(2,l),'FILTERED');
end
M;
N;
eNeuSnips=length(unique(M));
SnipSnips=length(unique(N));
if(eNeuSnips>SnipSnips)
    N=M;
end


catches = [];
Data = [freq(1,1:y); N]; Data2 = Data; % aggregates important data
for i = y:-1:1 % deletes catch trials in response map
    if levl(1,i) == 0
        catches = [catches Data2(2,i)]; % stores catch trial spike counts to establishing baseline rate
        Data2(:,i) = [];
    end
end
OutData1=Data2(1,:); 
OutData2=Data2(2,:);
a=length(OutData1);    
if length(OutData1)<101
    for j=a:101
        OutData1(j)=NaN;
        OutData2(j)=NaN;
    end
end

% baseline = mean(catches) % calculates baseline spontaneous firing rate

TT.SetGlobalV('RespectOffsetEpoc',0);
TT.SetEpocTimeFilterV('Freq',-.2,0);
trs = TT.GetValidTimeRangesV;
[x,y] = size(levl);
N2 = zeros(1,y);
if y > 110
    y = 110;
end
[q,p] = size(trs);
if y > p
    y=p;
end
N2 = zeros(1,y);
for l = 1:y    
    N2(l) = TT.ReadEventsV(100000,'Snip',channel,0,trs(1,l),trs(2,l),'FILTERED');
end

baseline=mean(N2);
TT.SetGlobalV('RespectOffsetEpoc',1);


TT.CloseTank
TT.ReleaseServer
close(figurea);


% figure
% plot(Data2(1,:),Data2(2,:),'k*')
% hold on
% plot(BFguess,baseline,'or')
% set(gca,'xscale','log');
% set(gca,'xlim',[min(freq(1,:)) max(freq(1,:))])
end
