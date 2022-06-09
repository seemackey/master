%  
% 
 function [] = MonkeyDiscPlotter(Monkey, Date, Num)
%% Create Handle; Import and Process Data
Monkey='Isildur';
Date='200831i';
Num='1';
Num=num2str(Num);
Date1=Date(1:6);


MyTank = ['Z:\' Monkey ' data\Behavior\ex' Date1 '\tank' Date];
MyBlock = ['\OurData-' Num];

data = TDTbin2mat(MyTank MyBlock)

[freq, freq2, DeltaStim, Responses] = DiscDataExtraction(MyTank,MyBlock)
% 
%%%%%%%%%%%%%%%%%%% Frequency plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check =  length(unique(freq(1,:))) 
% 
% if  check < 3
% 
% freqOut(1,:)= (freq(1,:));
% freqOut(2,:)=(freq2(1,:));
% freqOut(3,:)=Responses(1,:);
% 
% else
%     
% freqOut(2,:)= (freq(1,:));
% freqOut(1,:)=(freq2(1,:));
% freqOut(3,:)=Responses(1,:);
% 
% end
% 
% for i = 1:length(freqOut(1,:))
% diff(1,i) = (freqOut(2,i) - freqOut(1,i)) / freqOut(1,i);
% end
% 
% uniqueDiff = unique(diff);
% 
% for j = 1:length(uniqueDiff(1,:))
%     freqDiff(1,j) = uniqueDiff(1,j);
%     X = find(diff ==(uniqueDiff(1,j)));
%     l = 1;
%     for k = 1:length(X(1,:))
%      
%      rr(1,l) = Responses(1,X(1,k))
%      l = l+1;
%     end
%     cor = find(rr == 1);
%     TF = isempty(cor)
%     if TF == 1
%          ProbCor(1,j)=0
%     else
%     ProbCor(1,j) = length(cor(1,:)) / length(rr(1,:));
%     end
% end 
%   
% 
% figure
% semilogx(freqDiff,ProbCor)
% xlabel('delFreq')
% clc
% if check < 2
% display(['Frequecy1 is constant: ',num2str(freqOut(1,1))])
% else
% display(['Frequecy2 is constant: ',num2str(freqOut(1,1))])
% end
% 
% 
% 
% display(['Released for ', num2str(ProbCor(1,1)*100),' % of the catch trials'])
% display(['Trials Completed: ',  num2str(length(freqOut(1,:)))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot  tone levels %%%%%%%%%%
% check =  length(unique(freq(1,:))) 
% V0 = 0.79/(10^(74/20));   %volatge room 40
% if  check < 3
% 
% ToneOut(1,:)= (freq(1,:));
% ToneOut(2,:)=(freq2(1,:));
% ToneOut(3,:)=Responses(1,:);
% 
% else
%     
% ToneOut(2,:)= (freq(1,:));
% ToneOut(1,:)=(freq2(1,:));
% ToneOut(3,:)=Responses(1,:);
% 
% end
% 
% ToneOut(1,:) = 20*log10((ToneOut(1,:))/V0);
% ToneOut(2,:) = 20*log10((ToneOut(2,:))/V0);
% 
% for i = 1:length(ToneOut(1,:))
% diff(1,i) = (ToneOut(2,i) - ToneOut(1,i))         %/ freqOut(1,i);
% end
% 
% uniqueDiff = unique(diff);
% 
% for j = 1:length(uniqueDiff(1,:))
%     ToneDiff(1,j) = uniqueDiff(1,j);
%     X = find(diff ==(uniqueDiff(1,j)));
%     l = 1;
%     for k = 1:length(X(1,:))
%      
%      rr(1,l) = Responses(1,X(1,k))
%      l = l+1;
%     end
%     cor = find(rr == 1);
%     TF = isempty(cor)
%     if TF == 1
%          ProbCor(1,j)=0
%     else
%     ProbCor(1,j) = length(cor(1,:)) / length(rr(1,:));
%     end
% end 
%   
% 
% figure
% semilogx(ToneDiff,ProbCor)
% xlabel('delTone dB')
% clc
% if check < 2
% display(['Tone1 is constant: ',num2str(ToneOut(1,1))])
% else
% display(['Tone2 is constant: ',num2str(ToneOut(1,1))])
% end
% 
% 
% 
% display(['Released for ', num2str(ProbCor(1,1)*100),' % of the catch trials'])
% display(['Trials Completed: ',  num2str(length(ToneOut(1,:)))]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot  Modulation Depth %%%%%%%%%%
% V0 = 0.79/(10^(74/20));   %volatge room 40
check =  length(unique(freq(1,:))) 

if  check < 2

MDOut(1,:)= (freq(1,:));
MDOut(2,:)=(freq2(1,:));
MDOut(3,:)=Responses(1,:);

else
    
MDOut(2,:)= (freq(1,:));
MDOut(1,:)=(freq2(1,:));
MDOut(3,:)=Responses(1,:);

end

for i = 1:length(MDOut(1,:))
diff(1,i) = (MDOut(2,i) - MDOut(1,i))         %/ freqOut(1,i);
end

uniqueDiff = unique(diff);

for j = 1:length(uniqueDiff(1,:))
    ToneDiff(1,j) = uniqueDiff(1,j);
    X = find(diff ==(uniqueDiff(1,j)));
    l = 1;
    for k = 1:length(X(1,:))
     
     rr(1,l) = Responses(1,X(1,k))
     l = l+1;
    end
    cor = find(rr == 1);
    TF = isempty(cor)
    if TF == 1
         ProbCor(1,j)=0
    else
    ProbCor(1,j) = length(cor(1,:)) / length(rr(1,:));
    end
end 
  

figure
plot(ToneDiff,ProbCor)
xlabel('delMD')
clc
if check < 2
display(['MD1 is constant: ',num2str(MDOut(1,1))])
else
display(['MD2 is constant: ',num2str(MDOut(1,1))])
end



display(['Released for ', num2str(ProbCor(1,1)*100),' % of the catch trials'])
display(['Trials Completed: ',  num2str(length(MDOut(1,:)))]);




% End Session
%     Close tank and end session
%     invoke(TT, 'CloseTank');
%     invoke(TT, 'ReleaseServer');
%     close(f1); % closes handle function created in DataExtraction
%     Frequency=(mode(freq(1,:)))
%     display (Frequency)
%     display (threshold)
%     %%%Suppress the following line if not saving data
%   OutputImportantData(Monkey,Date,Responses,freq,freq2);
%     %%%%
end

