
function [] = MonkeyCurrentTankPlotter(Tank,Block)


 MyBlock = ['~OurData-' Block];
 MyTank=Tank;

[freq levl nlvl Responses TT f1] = DataExtraction(MyTank,MyBlock); % pulls out data

[TimeOfResponse TOR Responses guesscount] = ElimDuplTrials(Responses, levl); % Eliminate Duplicate Trials

RecordFreq = freq(1,1);
NLevel = nlvl(1,1);

[freq, levl, nlvl, Responses, TimeOfResponse] = EndTrialData(freq, levl, nlvl, Responses, TimeOfResponse); % Deletes trials after data ended



[ToneLevel alpha] = TLevlFindNOrder(levl); % Pulls out Tone Levels and orders them

[ProbabilityCorrect LevelCorr] = ResponsesNHits(levl,alpha,Responses,ToneLevel); % Pull out Responses to tone levels

%% Convert sound levels to dB SPL
LevelCorr(1,:);

alpha = length(LevelCorr(1,:));

LevelsUsed = dBtoSPL(ToneLevel, alpha); %converts tone levels from voltages to dbSPL

ProbabilityCorrect(2,:) = LevelsUsed;

threshold = ThreshCalc(ProbabilityCorrect); %calculates psychometric threshold
%% plot
figure
title('Probability Correct')
plot(ProbabilityCorrect(2,2:end),ProbabilityCorrect(1,2:end),'k-')
%  figure
%  title('Percent Correct')
%  plot(ProbabilityCorrect(2,2:end),LevelCorr(5,2:end),'r-')

%% End Session
% Close tank and end session
invoke(TT, 'CloseTank');
invoke(TT, 'ReleaseServer');
close(f1); % closes handle function created in DataExtraction
Frequency=(mode(freq(1,:)));
display (Frequency)
display (threshold)


end

