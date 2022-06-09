function MonkeysDoLotsOfStuff  %creating this function to do a LOT of analysis at the touch of a button... and a preset list of things to do. 


Dirrectory='C:\Users\Kimaya\Desktop\ThresholdData\QuickAnalysis.txt';

fid=fopen(Dirrectory, 'a');
fclose all; 
fid=fopen(Dirrectory, 'r');
Everything=textscan(fid, '%s%s%f%s');
fclose all;

iMonkeyName=Everything{1};
iMonkeyDataTank=Everything{2};
iTrialNum=Everything{3};
iRecordSite=Everything{4};
h = waitbar(0,'Please wait...');
for i=1:length(Everything{1})
    SaveThresh=1;
    
    
    MonkeyName=iMonkeyName{i};
    MonkeyDataTank=iMonkeyDataTank{i};
    TrialNum=iTrialNum(i);
    RecordSite=iRecordSite{i};
    MonkeyNandLevel(MonkeyName,MonkeyDataTank,TrialNum,RecordSite,SaveThresh);
    cab
    waitbar(i /length(Everything{1}))
end
close h

