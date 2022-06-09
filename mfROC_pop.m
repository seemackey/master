function [pc,pFA,pHit] = mfROC_pop(levelspikes,stdidx,bmf,VS)
% INPUT: "levelspikes" where each column is a trial, and each row is a condition, 
% and the value in
% each cell is number of spikes/or a spike timing value. 
% OUTPUTs:
% Hit rate, false alarm rate, pc (probability correct, AKA area under the
% ROC curve).
% 

ROCarray=levelspikes(stdidx:1:bmf,:); %subset of MFs we're going to analyze

%list of criterion values, change spacing depending on spike rate vs. VS
if VS==1
    critList = 0.01:0.05:max(levelspikes(:));
else
    critList = 1:1:max(levelspikes(:)); 
end



% calculate hits and false alarms at various criteria
% but preallocate first
pHit = zeros(length(critList),length(ROCarray(:,2)));
pFA =  zeros(length(critList),length(ROCarray(:,2)));

    for mfNum = 1:1:length(ROCarray(:,2))
        for critNum = 1:length(critList)
            criterion = critList(critNum);
            pHit(critNum,mfNum) = sum(ROCarray(mfNum,:)>criterion)/length(ROCarray(1,:));
            pFA(critNum,mfNum) = sum(ROCarray(1,:)>criterion)/length(ROCarray(1,:));
        end
    end

% plot ROC (HRs and FAs from last for loop)
plot(pFA,pHit,'-x');

% loop through and get the area under the ROC curve ("pc")
mfNum=1;
    for mfNum = 1:1:length(ROCarray(:,2))
        pc(mfNum,1) = -trapz(pFA(:,mfNum),pHit(:,mfNum));
    end

end