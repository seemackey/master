
function [ProbabilityCorrect,LevelCorr] = MonkeyTankPlotter2(freq, levl, nlvl, Responses,TT)


[TimeOfResponse TOR Responses guesscount] = ElimDuplTrials(Responses, levl); % Eliminate Duplicate Trials

RecordFreq = freq(1,1);
NLevel = nlvl(1,1);

[freq, levl, nlvl, Responses, TimeOfResponse] = EndTrialData(freq, levl, nlvl, Responses, TimeOfResponse); % Deletes trials after data ended

[ToneLevel alpha] = TLevlFindNOrder(levl); % Pulls out Tone Levels and orders them

[ProbabilityCorrect LevelCorr] = ResponsesNHits(levl,alpha,Responses,ToneLevel); % Pull out Responses to tone levels

%% Convert sound levels to dB SPL
alpha = length(LevelCorr(1,:));

LevelsUsed = dBtoSPL(ToneLevel, alpha); %converts tone levels from voltages to dbSPL

ProbabilityCorrect(2,:) = LevelsUsed;



end

