function [ProbabilityCorrect LevelCorr] = ResponsesNHits2(levl,alpha,Responses,ToneLevel)

% Pull out Responses to tone levels
% False alarm rate (wrong at 0 level)
% Hit rate (Right at a given level)
AllLevels=levl(1,:);
LevelCorr = zeros(5,alpha); %creates matrix to store responses
LevelCorr(1,:)=ToneLevel;
for i=1:alpha %nested loop to look at put all responses with correct tone level
    for j=1:length(AllLevels)
        if AllLevels(j) ~= 0
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
end


for j = 1:length(AllLevels)
    if AllLevels(j) == ToneLevel(1)
        LevelCorr(4,1)=LevelCorr(4,1)+1; % total trials at that level    
        if Responses(1,j) == 2
            LevelCorr(3,1) = LevelCorr(3,1) + 1;
        elseif Responses(1,j) == 3
            LevelCorr(2,1) = LevelCorr(2,1) + 1;
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

% for i=1:length(LevelCorrLow) %if there are 5 or less trials at a given tone level they are all deleted
%     LevelCorr(:,LevelCorrLow(i))=[];
%     ToneLevel(LevelCorrLow(i)) = [];
% end

alpha = length(ToneLevel);
% Calculate p(c) = ((hits + correct rejectios) / total trials
ProbabilityCorrect = zeros(1,alpha);
for i = 2:alpha
    ProbabilityCorrect(1,i) = 0.5*(LevelCorr(5,i)-LevelCorr(5,1))+0.5;
end

end