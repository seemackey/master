% Written by Andy Hrnicek under guidance of P.I. Ram Ramachandran, PhD.
% Fall 2010
% 
% this function calculates hit rate and choice probability of the function.
% It pulls out and bins the data according to ToneLevel of the trial and
% the response (1 or 2) and then calculates false alarm rate and hit rate,
% from there it is converted into choice probability.
%
%%% INPUTS
% levl - tone level by trial
% alpha - number of tone levels used in this block
% Responses - response by trial (1 or 2)
% ToneLevel - tone levels in ascending order in volts
%%% OUTPUTS
% ProbabilityCorrect - choice probability values in ascending order of tone
% level
% LevelCorr - Aggregate data of responses by tone level. Row1 - Tone Level
% in volts, Row2 - number of incorrect responses, Row3, number of correct
% responses, row4 - total trials at that tone level, row5 - hit rate
% ToneLevel - tone levels in ascending order in volts

function [ProbabilityCorrect LevelCorr ToneLevel] = ResponsesNHits(levl,alpha,Responses,ToneLevel)
% ProbabilityCorrect(1,:) - [FA, choice probability values for tone levels]
% ProbabilityCorrect(2,:) - [Tone levels in dBSPL(assigned in TankPlotter)]
%
% LevelCorr(1,:) - ToneLevels in volts
% LevelCorr(2,:) - count of Response = 2 for that tone level
% LevelCorr(3,:) - count of Response = 1 for that tone level
% LevelCorr(4,:) - total trials for that tone level
% LevelCorr(5,:) - Hit Rates
% 
% ToneLevel - array of tone levels in voltages

% Pull out Responses to tone levels
% False alarm rate (wrong at 0 level)
% Hit rate (Right at a given level)



AllLevels=levl(1,:);
LevelCorr = zeros(5,alpha); %creates matrix to store responses
LevelCorr(1,:)=ToneLevel;
for i=1:alpha %nested loop to look at put all responses with correct tone level
    for j=1:length(Responses)
        if ToneLevel(i)==AllLevels(j)
            LevelCorr(4,i)=LevelCorr(4,i)+1; % total trials at that level
            if Responses(1,j)==1
                LevelCorr(3,i)=LevelCorr(3,i)+1; %Row 3 is Corr=1
            elseif Responses(1,j)==2
                LevelCorr(2,i)=LevelCorr(2,i)+1; %Row 2 is Corr=2
            end
        end
    end
end

% Calculate hit rate for each tone level
LevelCorrLow = [];
for i=1:alpha
    LevelCorr(5,i) = LevelCorr(3,i)/LevelCorr(4,i);
    if LevelCorr(4,i)<6
        LevelCorrLow = [i LevelCorrLow];
    end
end

for i=1:length(LevelCorrLow) %if there are 5 or less trials at a given tone level they are all deleted
    LevelCorr(:,LevelCorrLow(i))=[];
    ToneLevel(LevelCorrLow(i)) = [];
end

alpha = length(ToneLevel);
% Calculate: p(c) = .5*(HR(tlevl)-FA)+.5
ProbabilityCorrect = zeros(1,alpha);


for i=1:alpha
    if LevelCorr(5,i)==1
        LevelCorr(5,i)=.99;
    elseif LevelCorr(5,i)==0
        LevelCorr(5,i)=.01;
    end
end


PercentCorrect=vertcat(LevelCorr(1,:),LevelCorr(5,:));
dprime1=norminv(LevelCorr(5,:))-norminv((LevelCorr(5,1)));
ProbabilityCorrect=dprime1; % substituting dprme for PC temporarily for duration paper
ProbabilityCorrect=normcdf(dprime1./2); % convert to prob correct
ProbabilityCorrect(1)=NaN;

% for i = 2:alpha    
%     ProbabilityCorrect(1,i) = .5*(LevelCorr(5,i)-LevelCorr(5,1))+.5; % the other equation option is less reliable because we don't do the same number of trials for all tone levels
% end
end