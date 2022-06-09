

function [] = MonkeyTankPlotter(Monkey, Date, Num, x, Type)
%% Create Handle; Import and Process Data



Num=num2str(Num);
Date1=Date(1:6);
% sets the location of the block and tank to be accessed, this is different
% on all computers and needs to be changed for it to work properly.
    if Type==1
        MyTank = ['G:\' Monkey ' data\Behavior\ex' Date1 '\tank' Date];
    else
        MyTank = ['G:\' Monkey ' data\Neurobehavior\' x '\ex' Date1 '\tank' Date];
    end
    MyBlock = ['~OurData-' Num];
    

    
    
    
[freq,  levl, nlvl, nLvls,  Responses, TT, f1] = DataExtraction(MyTank,MyBlock); % pulls out data

% [TimeOfResponse TOR Responses guesscount] = ElimDuplTrials(Responses, levl); % Eliminate Duplicate Trials
% Responses
%Row 1 = correct or incorrect
%Row 2 = time of trial ending or corr
% levl
%Row 1= Voltage level of tone
%Row 2= onset time of levl

% freq
%Row 1= Frequency
%Row 2= Onset time of Freq
%Row 3= Offset time of Freq




% Use this code to help plot attention experiment

% freqOut = zeros(2,size(levl,2));
freqOut(1,:)=round(freq(1,:));
freqOut(2,:)=Responses(1,:);

Uniq=(unique(freqOut(1,:)));
Uniq(2:3,:)=0;

if length(Uniq(1,:))>3 %if more than one frequency, it will do this loop (attention)
    NewLength=floor(length(freqOut)/2);
    for(i=1:NewLength)
        
        Probes(:,i)=freqOut(:,2*i);
        
    
    end
    Probes
    length(unique(Probes))
    if length(unique(Probes(1,:)))==1
        for(i=1:NewLength)
            
            Probes(:,i)=freqOut(:,2*i-1);
            
            
        end
    end
    CFreq=mode(freqOut(1,:));   %find the center frequency of attention
    for i=1:length(Uniq(1,:))   
        
        i;
        a=find(Probes(1,:)==Uniq(1,i)); %go through by each unique frequency, and pool data
        b=find(Probes(2,a)==1);
        c=find(Probes(2,a)==2);
        d=length(b)+length(c);
        Uniq(2,i)=length(b);
        Uniq(3,i)=length(b)/d;
    end
    

    semilogx(Uniq(1,:),Uniq(3,:))


else

    
    
    
    RecordFreq = freq(1,1);
    NLevel = nlvl(1,1);
    
    
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
   
    threshold = ThreshCalc(ProbabilityCorrect); %calculates psychometric threshold
%     LevelsUsed
    LevelCorr(5,:)
    
    
    AllArray(1,:)=freq(1,:); % tone
    AllArray(2,:)=round(dBtoSPL(levl(1,:),length(levl(1,:))));%voltage
    %AllArray(2,:)=levl(1,:);
    AllArray(3,:)=freq(2,:); %time at which the tone starts playing
    AllArray(4,:)=Responses(1,:);%correct resp
    AllArray(5,:)=Responses(2,:)-freq(2,:);%RT
    AllArray(6,:)=nlvl (1,:);%Noise Amp
   
    
    AllArray;
    
    Array=(AllArray');
    AmpArray=nLvls;
save Array

   
    
   
   
    %try    RwdVsUnrwd(AllArray,threshold);
    %catch
    %end
    
    %plot
    
%     AllArray(1,:)=freq(1,:);
%     AllArray(2,:)=dBtoSPL(levl(1,:),length(levl(1,:)));
%     AllArray(3,:)=freq(2,:);
%     AllArray(4,:)=Responses(1,:);
%     AllArray(5,:)=Responses(2,:)-freq(2,:);
%     AllArray;
%     
%     DidntRespondInd=find(AllArray(5,:)>(max(AllArray(5,:))-.0001));
%     
%     Unique_dB=unique(AllArray(2,:));
%     NumberUnique_dB=length(Unique_dB);
% 
% 
%     
%     
%     
%     
%     j=0;
%     
%     if ~isnan(threshold)
%         First=find(Unique_dB>threshold); %find first occurence above threshold
%         for i=First(1):First(end)
%             j=j+1;
%             Indeces=find(AllArray(2,:)==Unique_dB(i));
%             CorrectIndeces=find(AllArray(4,Indeces)==1);
%             CorrectIndeces=Indeces(CorrectIndeces);
%             MeanRxnTime(1,j)=Unique_dB(i);
%             MeanRxnTime(2,j)=mean(AllArray(5,CorrectIndeces));
%             AllArray(5,CorrectIndeces)=(AllArray(5,CorrectIndeces)-mean(AllArray(5,CorrectIndeces)))/std(AllArray(5,CorrectIndeces));
%         end
%         %Now have adjusted rxn times for all correct trials over threshold
%         
%         First=find(Unique_dB<=threshold); %find first occurence above threshold
%         j=0;
%         for i=First(1):First(end)
%             j=j+1;
%             Indeces=find(AllArray(2,:)==Unique_dB(i));
%             CorrectIndeces=find(AllArray(4,Indeces)==1);
%             CorrectIndeces=Indeces(CorrectIndeces);
%             MeanLRxnTime(1,j)=Unique_dB(i);
%             MeanLRxnTime(2,j)=mean(AllArray(5,CorrectIndeces));
%             AllArray(5,CorrectIndeces)=(AllArray(5,CorrectIndeces)-mean(AllArray(5,CorrectIndeces)))/std(AllArray(5,CorrectIndeces));
%         end

    
    figure
    title('Probability Correct')
    plot(ProbabilityCorrect(2,2:end),ProbabilityCorrect(1,2:end),'k-')
    % figure
    % title('Percent Correct')
    % plot(ProbabilityCorrect(2,2:end),LevelCorr(5,2:end),'r-')
    lvl(:,1) = ProbabilityCorrect(2,2:end) 
    %HR(:,2) = LevelCorr(5,2:end) comment this in if you want hit rate
    %HR = flip(HR)
    %FA=ProbabilityCorrect(1,:)
    PC=ProbabilityCorrect(1,2:end);
    PC = PC'
    
    
    %% End Session
    % Close tank and end session
    invoke(TT, 'CloseTank');
    invoke(TT, 'ReleaseServer');
    close(f1); % closes handle function created in DataExtraction
    Frequency=(mode(freq(1,:)));
    display (Frequency)
    display (threshold)
    %%%%Suppress the following line if not saving data
  OutputImportantData(LevelCorr,ProbabilityCorrect,Monkey,Date,Num,Responses,levl,freq,threshold);
    %%%%
end
end
