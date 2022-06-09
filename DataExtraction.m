% Data Extraction
% Written by Ram Ramachandran, adapted to current form by Andy Hrnicek
%
% This code accesses the TDT data and stores it in variables matlab can
% use.
%
%%% INPUTS
% MyTank - string of pathway to data tank of interest
% MyBlock - string describing block of interest
%%% OUTPUTS
% freq - trial frequency and time stamps
% levl - trial tone level and time stamps
% nlvl - trial noise level and time stamps
% Responses - response by trial (1 or 2)
% TT - TTank.x figure handle
% figure1 - figure handle for TTank.x handle

function [freq levl nlvl nLvls Responses TT figure1] = DataExtraction(MyTank,MyBlock)
% this code creates the data handle and accesses data in the tank and
% stores it in the variables 'freq' - tone frequencies, 'levl' - tone level
% played, 'nlvl' - noise level played, and 'Responses' - 1 or 2 for correct
% and incorrect response by the animal

%Open tank/block for reading
figure1 = figure;
TT = actxcontrol('TTank.X');
invoke(TT,'ConnectServer','Local','Me');
invoke(TT,'OpenTank',MyTank,'R');
invoke(TT,'SelectBlock',MyBlock);
invoke(TT,'ResetFilters');

Responses = invoke(TT, 'GetEpocsV', 'Corr', 0, 0, 10000); %puts the 'Corr' variable into an x,1 matrix with lots of zeros
Responses(1,:); 

% Specify channel
% code = invoke(TT, 'GetEventCodes', 0);
% start = invoke(TT,'CurBlockStartTime')
% formstart = invoke(TT,'FancyTime',start ,'D/O/Y H:M:S.U W')
% notes = invoke(TT,'CurBlockNotes')
channel = 1;
invoke(TT,'SetGlobalV','Channel',channel);
freq = invoke(TT, 'GetEpocsV', 'Freq', 0, 0, 10000);

levl = invoke(TT, 'GetEpocsV', 'Levl', 0, 0, 10000);

nLvls= invoke(TT, 'GetEpocsV', 'NseL',0, 0, 30000); 
nlvl = invoke(TT, 'GetEpocsV', 'nlvl', 0, 0, 10000);

if length(Responses(1,:))<length(freq(1,:))
  
   
    freq=freq(:,1:end-1);

    levl=levl(:,1:end-1) ;
    nlvl=nlvl(:,1:end-1);
    nLvls=nLvls(:,1:end-1);
end


% invoke(TT,'SetFilterWithDescEx',['Freq=' num2str(f1)]);
% wave = invoke(TT, 'ReadWavesV', 'Wave');
% pdec = invoke(TT, 'ReadWavesV', 'PDec');

% wave = TT.ReadEventsV(100000, 'Wave', channel, 0, 0, 0, 'ALL');
% pdec = TT.ReadEventsV(100000, 'PDec', channel, 0, 0, 0, 'ALL');
% 
% wave = TT.ParseEvV(0,wave);
% pdec = TT.ParseEvV(100000, 'PDec', channel, 0, 0, 0, 'ALL');

% TT.CreateEpocIndexing
% N = TT.ReadEventsV(1000000, 'Snip', channel, 0, 0.0, 0.0, 'FILTERED');
% TS = TT.ParseEvInfoV(0, N, 6);
% F = TT.QryEpocAtV('Levl', TS(1075), 0);
% F1 = TT.QryEpocAtV('Levl', TS(1075), 1);
% F2 = TT.QryEpocAtV('Levl', TS(1075), 2);
% F3 = TT.QryEpocAtV('Levl', TS(1075), 3);

% error checking to make sure that there is adequate data in the block for
% analysis and to make sure that the block exists in the location specified
MyBlock

if strcmp(num2str(freq(1,1)),'NaN')
    MyDate = MyTank((end-5):end);
    fprintf('Block %s, Date %s\r',MyBlock,MyTank);
    error('Tank_Location:NotFound','Data Tank may not exist or specified location is incorrect.')
end
[m,n] = size(freq);
if n < 50
    MyDate = MyTank((end-5):end);
    fprintf('Block %s, Date %s\r',MyBlock,MyDate);
    error('Short_Tank:Short','Data Tank may be incomplete or a short trial.')
end 

end