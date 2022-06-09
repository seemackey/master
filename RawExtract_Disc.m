function [All] = RawExtract(path,monkey)

%% This function is good for looping through so you can get all trials in an
%  array and combine them for analysis or modeling

%% I/O
% INPUTS %
% path - the path to the data
% monkey - self explanatory

% OUTPUTS %
% All - an array with all data by trial

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
    % delete last trial to compensate
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
    
   
end