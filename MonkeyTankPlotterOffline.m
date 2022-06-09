close all;
clear all;

 
 %% ENTER PATH INFO HERE
Monkey='Opal';
Date='220503o';
Num='2';
Num=num2str(Num);
Date1=Date(1:6);
lvlup=0;

path = ['/Volumes/ramlab2018/' Monkey ' data/Behavior/ex' Date1 '/tank' Date '/OurData-' Num];

%path = 'Y:\Alpha data\Behavior\ex170116\ex170116\tank170116a\OurData-6'
%path = 'T:\Gandalf data\Behavior\ex191216\tank191216g\OurData-4'; %syn 040
%path = 'Y:\Alpha data\Behavior\ex181026\tank181026a\OurData-8'; %ramlab
%path = 'T:\Haldir data\Behavior\ex191218\tank191218h\OurData-1'; %haldir teba
% path = 'T:\Bilbo data\Behavior\ex190618\tank190618b\OurData-5'; %ram2018
%path = 'T:\Aragorn data\Behavior\ex191209\tank191209a\OurData-6';
%path = 'Y:\Bravo data\Behavior\ex190617\tank190617b\OurData-4'; %ramlab

% threshold plus these values for rt allocation in outputimptdata function
min=20;
max=26;


Data = TDTbin2mat(path);

% Check if responses are equal to presentations. (Caught blocks with one less response)

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
    All(9,:) = Data.epocs.Freq.offset';
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
    All(9,:) = Data.epocs.Freq.offset';
end

% get duration of stimulus
Durs=Data.epocs.Freq.offset-Data.epocs.Freq.onset;
Duration=median(Durs); % only use this as your dur if the dur is constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FALSE ALARM RATE CALC. 
% counts up # of catch trials
ct=0;
for b=1:length(All)
    if All(3,b)==0
        ct=ct+1;
    else
        continue
    end
end
% % % counts up # of false alarms
fa=0;
for bb=1:length(All)
    if All(4,bb)==1 && All(3,bb)==0
        fa=fa+1;
    else 
        continue
    end
end
% % % divide them both to get false alarm rate
FArate=fa/ct

freq = All(1,:);
freq(2,:) = All(8,:);
freq(3,:) = All(9,:);
levl =  All(3,:);
levl(2,:) = All(6,:);
nlvl = All(2,:);
nlvl(2,:) = All(7,:);
Responses = All(4,:);
Responses(2,:) = All(5,:);

if length(Responses(1,:))<length(freq(1,:))
  
   
    freq=freq(:,1:end-1);

    levl=levl(:,1:end-1) ;
    nlvl=nlvl(:,1:end-1);
    nLvls=nLvls(:,1:end-1);
end
    
    

    RecordFreq = freq(1,1);
    NLevel = nlvl(1,1);
%     AmpN=nLvls(1,1);
    
%     [freq, levl, nlvl, Responses, TimeOfResponse] = EndTrialData(freq, levl, nlvl, Responses, TimeOfResponse); % Deletes trials after data ended
    
    
    [ToneLevel alpha] = TLevlFindNOrder(levl); % Pulls out Tone Levels and orders them
    ToneLevel;
    alpha;
    
    [ProbabilityCorrect LevelCorr] = ResponsesNHits(levl,alpha,Responses,ToneLevel); % Pull out Responses to tone levels
    
    %% Convert sound levels to dB SPL
    LevelCorr(1,:);
    
    alpha = length(LevelCorr(1,:));
    
    LevelsUsed = dBtoSPL(ToneLevel, alpha); %converts tone levels from voltages to dbSPL
    
    ProbabilityCorrect(2,:) = LevelsUsed;
    
    threshold = ThreshCalc(ProbabilityCorrect) %calculates psychometric threshold
   disp(LevelCorr); %levels
    
   %% psych function/accuracy plotting, printing etc.
   figure
    title('Probability Correct')
    plot(ProbabilityCorrect(2,2:end),ProbabilityCorrect(1,2:end),'Marker','o')
    %
    
    
    % figure
    % title('Percent Correct')
    %plot(ProbabilityCorrect(2,2:end),LevelCorr(5,2:end),'r-')
    
    %FArate=LevelCorr(5,1) %FALSE ALARM RATE %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    HR(:,1) = ProbabilityCorrect(2,2:end);
    HR(:,2) = LevelCorr(5,2:end);
    PC = ProbabilityCorrect(1,2:end);
    PC = PC'; % probability correct
    %% End Session
    % Close tank and end session
    
    Frequency=(mode(freq(1,:)));
     display (Frequency)
    %display (threshold)
    %%%%Suppress the following line if not saving data
    
    % setting levels for RT allocation
    RTmin=threshold+min;
    RTmax=threshold+max;
  OutputImportantData(LevelCorr,ProbabilityCorrect,Responses,levl,freq,threshold,RTmin,RTmax);
  
  
  x = LevelsUsed(1,2:end)';
  y = ProbabilityCorrect(1,2:end)';
    %%%%
    
 %%   %%%%%%%%%%%%%%%%%Creates Weibull Fit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    createFitWB(HR,PC,ProbabilityCorrect(2,2),lvlup)
    
    % saves the psych function figure to current directory
    saveas(figure(3),[pwd Monkey Date Num '.fig']);
    
    
    %thresh=threshold-40;
   
    
    
    
    %RTvCount(levl,Responses) %what does this do?
 %lvls=HR(:,1); % these are the actual levels used, translated back down
 

