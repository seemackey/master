function [All] = RawExtract(path,monkey)

%% This function is good for looping through so you can get all trials in an
%  array and combine them for analysis or modeling

%% I/O
% INPUTS %
% path - the path to the data
% monkey - self explanatory

% OUTPUTS %
% All - an array with all data by trial

%% extract data struct
Data = TDTbin2mat(path);

%% Check if responses are equal to presentations
LF =  length(Data.epocs.Freq.data(:,1));
LC = length(Data.epocs.Corr.data(:,1));
    if isequal(LF,LC)
        % ALL in one matrix
        All(1,:) = Data.epocs.Freq.data';
        All(2,:) = Data.epocs.nlvl.data';
        All(3,:) = Data.epocs.Levl.data';
        All(4,:) = Data.epocs.Corr.data';
        All(5,:) = Data.epocs.Corr.onset';
        All(6,:) = Data.epocs.Levl.onset';
        All(7,:) = Data.epocs.nlvl.onset';
        All(8,:) = Data.epocs.Freq.onset';
    else

        %Delete last trial to compensate one less response
        Data.epocs.Freq.data(LF,:) = [];
        Data.epocs.nlvl.data(LF,:) = [];
        Data.epocs.Levl.data(LF,:) = [];


        % ALL in one matrix
        All(1,:) = Data.epocs.Freq.data';
        All(2,:) = Data.epocs.nlvl.data';
        All(3,:) = Data.epocs.Levl.data';
        All(4,:) = Data.epocs.Corr.data';
        All(5,:) = Data.epocs.Corr.onset';
        All(6,:) = Data.epocs.Levl.onset';
        All(7,:) = Data.epocs.nlvl.onset';
        All(8,:) = Data.epocs.Freq.onset';
    end
    %(:,3)=dBtoSPL(All(:,3),length(All(:,3))); convert to dB?
    All(5,:)=-(All(6,:)-All(5,:))*1000; % subtract stim onset and RT, and convert to ms
    All(7,:)=monkey;
    
     % get duration of stimulus
    Durs=Data.epocs.Freq.offset-Data.epocs.Freq.onset;
    Duration=median(Durs); % only use this as your dur if the dur is constant
    All(8,:)=Duration;
    
   
end