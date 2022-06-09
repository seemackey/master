% loop through the samsmodlong pooling fxn to get a count of how many times
% threshold is met by the populations of neurons
function [ChoiceProbs,Grand_CP,p_value,grandcpboot] = CP_Calc (data, dur)
format long g; 


 %% Loop through data and run ROC/CP Analysis on each of the neurons 

     %Initialize all variables to input into CP function to match format in MonkeyNandLevel loop 
 
        tlvls = (unique(data.epocs.Levl.data))';
        channel = 1;
        
        Responses = [
            (data.epocs.Corr.data)';
            (data.epocs.Corr.onset)';
                ];
          
        freq = [
            (data.epocs.Freq.data)';
            (data.epocs.Freq.onset)';
            (data.epocs.Freq.offset)';
                ];
        global freqout
        freqout=rot90(freq(1,:));
        levl = [
            (data.epocs.Levl.data)';
            (data.epocs.Levl.onset)';
                ];
        nlvl = [
             zeros(1,length(data.epocs.Levl.onset));
            (data.epocs.Levl.onset)';  
                ];   

            part1=length(freq(1,:));
            part2=length(unique(freq(1,:)));
            BFused=mode(freq(1,:));

            part4=min(freq(1,:));
            warn1str=num2str(part1);
            warn1str=['Only ',warn1str, ' Trials Recorded, Terminating Analysis'];
            if part1<70
                warndlg(warn1str, 'Warning!');
                disp(data.info.tankpath)
                close(figurea);
                levl=0;
                N=0;
                return
            elseif part2>3 && part4>0
                warndlg('Block selected is most likely Response Map', 'Warning!');
                disp(data.info.tankpath)
                close(figurea);
                levl=0;
                N=0;
                return
            end

            %tlvls = TLevlFindNOrder(levl)
            tlvls = unique(levl(1,:)); % contains tone levels used

            CheckIfIsMod=find(freq(1,:)==12321);
            if isempty(CheckIfIsMod)==0
                ModNoise=1;
                freq(1,CheckIfIsMod)=median(freq(1,:));
            else
                ModNoise=0;
            end

            %CHANGE NOISE LEVEL BASED ON THE ROOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ROCFreq=median(freq(1,:));  %the frequency for the ROC analysis
            ROCNoise=median(nlvl(1,:)); %ROC Noise Level
            if ROCNoise>0
                if room==40
                v0 = .79/(10^(74/20)); % On Hylebos-PC in room 040
                else
                v0=.79/(10^(82/20)); % On Yakima in Room 034
                end
                ROCNoise=round(20*log10((ROCNoise)/v0));
            elseif  ROCNoise==0
                ROCNoise=24;
            end

            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % psychometric function calc.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %ProbabilityCorrect=MonkeyTankPlotter2(freq,levl,nlvl,Responses,TT);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% importing trial info
            trs = [
                (data.epocs.Freq.onset)';
                (data.epocs.Freq.offset)'; 
                ]; %gets start and end time of each epoch 
            

            maxtime = 1/325;
            [~,y] = size(levl);
            [~,y1] = size(trs);
            if y > 400 % deletes extraneous trials
                y = 400;
                freq = freq(:,1:400);
                levl = levl(:,1:400);
                nlvl = nlvl(:,1:400);
            end

            if y>y1 % deletes extraneous trials
                y = y1;
                freq = freq(:,1:y1);
                levl = levl(:,1:y1);
                nlvl = nlvl(:,1:y1);
            end

            r = rem(y,10);
            if r>0 % deletes extraneous trialsclc
                y = y-r;
                freq = freq(:,1:y);
                levl = levl(:,1:y);
                nlvl = nlvl(:,1:y);
            end

            N = zeros(1,y);
            M=N;
            DeleteCount = zeros(1,y);
            SnipTimes = zeros(y,20);
            InstRates = zeros(y,20);

            %%
            %%%%%%%%%%%%%%%% CHANGE DURATION OF ANALYSIS WINDOW HERE %%%%%%%%%%%%%%
            WindowAnalysis = 1;
            start=0;
            room=40;
            % 
            if WindowAnalysis==1
                trs(1,:)=trs(1,:)+start; %changing start of trial 

                trs(2,:)=trs(2,:)-(0.2-(dur+start)); %changing end of trial
                windowsize=trs(2,1)-trs(1,1); %make sure the math is right
            else

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  extract the spike counts in between the onset and offset times for each tone
data = TDTfilter(data,'Levl','TIME',[0, 0.2]);
            
            if isfield(data.snips,'Snip')
                ts = [data.snips.Snip.ts];
            else
                ts = [data.snips.eNeu.ts];
            end
                counter = 0;

            for trs_ct = 1:y
                for spike_counter = 1:length(ts)
                    if (ts(spike_counter) >= trs(1,trs_ct)) && (ts(spike_counter) <= trs(2,trs_ct))
                        counter = counter +1;
                    end
                end
                N(trs_ct) = counter;
                counter = 0;
                %SnipTimes(trs_ct,1:N(trs_ct)) = 
            end

            %% spikes and spike times have been extracted, now sorting spike data


            Sorted = SortingByTime(levl, Responses);

            m=length(Sorted(1,:));
            for i=1:m-1
                Sorted(4,i)=N(i); % puts # spikes in each trial in the 4th row of "Sorted"
            end


            a=length(levl);
            for i=1:a
                levl(2,i)=N(i);
            end
            b=length(tlvls);
            ToneLevels=tlvls;
            LengthTlvls=b;



            %%%%%%%%%%%%%%%%%%
            %                   Build an array with Tone levels as the first column,
            %                   assuming that the tone is repeated 30 times max. If
            %                   this number is greater, make length the number of
            %                   repeats +1
            %%%%%%%%%%%%%%%%%%
            for i=2:35
                tlvls(i,:)=NaN;
            end



            %%%%%%%%%%%%%%%%%%
            %                   The following for loop sorts spike counts by tone from
            %                   the sorted array into two arrays, tlvlsCorrect which
            %                   corresponds to a correct pull from the monkey, and
            %                   tlvlsIncorrect which contains all spike data miss
            %                   trials
            %%%%%%%%%%%%%%%%%%


            tlvlsCorrect=tlvls;
            tlvlsIncorrect=tlvls;
            for i=1:b
                count1=2;
                count2=2;
                for j=1:a
                    if Sorted(1,j)==tlvls(1,i)
                        if Sorted(3,j)==1
                            tlvlsCorrect(count1,i)=Sorted(4,j);
                            count1=count1+1;
                        elseif Sorted(3,j)==2
                            tlvlsIncorrect(count2,i)=Sorted(4,j);
                            count2=count2+1;
                        end
                    end
                end
            end


            NANflag=0;
            for i=1:length(tlvlsCorrect(1,:))
                Check1=tlvlsCorrect(2,i);
                Check2=tlvlsIncorrect(2,i);
                if isnan(Check1)==1 && isnan(Check2)==1
                    NANflag=1;
                    RowIndex=i;
                end
            end
            if NANflag==1
                tlvlsCorrect(:,RowIndex)=[];
                tlvlsIncorrect(:,RowIndex)=[];
            end

            % sortin spike counts by tone level (in volts)
            for i=1:b
                c=2;
                for j=1:a
                    if tlvls(1,i)==levl(1,j)
                        tlvls(c,i)=levl(2,j);
                        c=c+1;
                    end
                end
            end
            if NANflag==1
                LengthTlvls=LengthTlvls-1;
                b=b-1;
                tlvls(:,RowIndex)=[];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%   Reshaping Tone level Data and Spitting Out Choice Probablity %%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%   Jackson Mayfield 08/30/21   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            [ChoiceProbs,Grand_CP,p_value,grandcpboot] = CP(tlvlsCorrect,tlvlsIncorrect,tlvls);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


end

function [LevelAndResponses] = SortingByTime(levl, Responses)

a=length(levl);
b=length(Responses);
j=1;

for i=1:a
    while(levl(2,i)>=Responses(2,j))
        j=j+1;
        if j>b
            LevelAndResponses=levl;
            return
        end
        
    end
    levl(3,i)=Responses(1,j);
    levl(5,i)=Responses(2,j)-levl(2,i);
end
LevelAndResponses=levl;
end
