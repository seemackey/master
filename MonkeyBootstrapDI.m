

function [WeightedDI]=MonkeyBootstrapDI(Correct,Incorrect)





%Need to make new array of correct and incorrect.... actually just need to
%:LKHF:LEKH:ALKHF:LKAHF 

% Dimension has number of correct in first row, and number of incorrect in second row. We can use this to select (C) and (I)
% AllDIs=zeros(1,20);
% for qwerty=1:20

a=length(Correct(1,:));
b=length(Correct(:,1));
Dimension=zeros(2,a);
DI=zeros(3,a);
DI(1,:)=Correct(1,:);
Greater=zeros(1,a);
Equal=zeros(1,a);

Correct(1,:)=-Correct(1,:);        %Making the first row of correct and incorrect 
Incorrect(1,:)=-Correct(1,:);      %inverted. 

SpikesOnly=sort(vertcat(Correct(2:end,:),Incorrect(2:end,:)));  %putting all spikes into one array, then sorting 

%Have two arrays, one of correct trials, for i=1:a %a is height correct/incorrect
for i=1:a
    for j=2:b %b is length correct/incorrect
        if isnan(Correct(j,i))==0
            Dimension(1,i)=Dimension(1,i)+1;
        end
        if  isnan(Incorrect(j,i))==0
            Dimension(2,i)=Dimension(2,i)+1;
        end
    end
end

%Now i have spikes in one array sorted by level, and number of correct and
%incorrect... just need randomize and put into "correct" and "incorrect"
%arrays again. 
NewCorrect=Correct*NaN;
NewCorrect(1,:)=Correct(1,:);
NewIncorrect=Incorrect*NaN;
NewIncorrect(1,:)=Incorrect(1,:);
NewDimension=Dimension*0;
%copied arrays and filled with NaNs
for x=1:length(Dimension(1,:))
    for y=1:sum(Dimension(:,x))
        Rando=ceil(sum(Dimension(:,x))*rand);
        if Rando>Dimension(1,x)
            NewDimension(2,x)=NewDimension(2,x)+1;
        else
            NewDimension(1,x)=NewDimension(1,x)+1;
        end
    end
end

% Now I SHOULD have a new bootstrapped Dimension array, should use this to
% make a new correct and incorrect array... 
SpikeAmounts=isnan(SpikesOnly);
SumSpikeAmounts=SpikeAmounts(1,:)*NaN;
for x=1:length(SpikeAmounts(1,:))
    SumSpikeAmounts(x)=find(SpikeAmounts(:,x),1,'first')-1; 
end
%SumSpikeAmounts lets me know how many real numbers are in the spike array
%for later resampling

for x=1:length(NewDimension(1,:))
    
    if NewDimension(1,x)>0
        for y=1:NewDimension(1,x)  %This loop is to resample Incorrect Trials, want to pull random data from SpikesOnly array
            
            NewCorrect(y+1,x)=SpikesOnly(ceil(rand*SumSpikeAmounts(x)),x);
        end
    end
    
    if NewDimension(2,x)>0
        for z=1:NewDimension(2,x) %This loop is to resample Correct Trials
            
            NewIncorrect(z+1,x)=SpikesOnly(ceil(rand*SumSpikeAmounts(x)),x);
        end
    end
end


%Now that i have a newly sorted incorrect and correct, lets START
%EVERYTHING over from the beginning, this time using NEW arrays... 
Correct=NewCorrect;
Incorrect=NewIncorrect;
a=length(Correct(1,:));
b=length(Correct(:,1));
DI=zeros(3,a);
DI(1,:)=Correct(1,:);
Greater=zeros(1,a);
Equal=zeros(1,a);
Dimension=zeros(2,a);



for i=1:a
    for j=2:b %b is length correct/incorrect
        if isnan(Correct(j,i))==0
            Dimension(1,i)=Dimension(1,i)+1;
        end
        if  isnan(Incorrect(j,i))==0
            Dimension(2,i)=Dimension(2,i)+1;
        end
    end
end


for i=1:a
    j=2;
    for l=2:b
        Comp1=Correct(j,i);
        for k=2:b
            if Comp1 > Incorrect(k,i)
                Greater(1,i)=1+Greater(1,i);
            end
            if Comp1 == Incorrect(k,i)
                Equal(1,i)=Equal(1,i)+1;
            end
        end
        j=j+1;
    end
end

SumWeight=0;
for i=1:a    
    DI(2,i)=(Greater(i)+.5*Equal(i))/(Dimension(1,i)*Dimension(2,i));
    if isnan(DI(2,i))==0;
        if i==1
            DI(3,i)=DI(2,i)*Dimension(1,i);
            SumWeight=Dimension(1,i)+SumWeight;
        else
            DI(3,i)=DI(2,i)*Dimension(2,i);
            SumWeight=Dimension(2,i)+SumWeight;
        end
    end
end

WeightedDI=sum(DI(3,:)/SumWeight);


% end
% 
% 
% AllDIs
% mean(AllDIs)
