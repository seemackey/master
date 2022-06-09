function [OnsetDelay] = MonkeyRespCummulative (Monkey, Date, Num, Loc)

% Monkey='Bravo';
% Date='111007b';
% Num='2';
% Loc='CN';
Num=num2str(Num);
Date1=Date(1:6);
MyTank = ['G:\' Monkey ' data\NeuroBehavior\' Loc '\ex' Date1 '\tank' Date];
MyBlock = ['~OurData-' Num];

figurea = figure;
TT = actxcontrol('TTank.X');
invoke(TT,'ConnectServer','Local','Me');
invoke(TT,'OpenTank',MyTank,'R');
invoke(TT,'SelectBlock',MyBlock);
invoke(TT,'ResetFilters');
TT.SetGlobalV('RespectOffsetEpoc',1);
Responses = invoke(TT, 'GetEpocsV', 'Corr', 0, 0, 10000); %puts the 'Corr' variable into an x,1 matrix with lots of zeros
channel = 1;
invoke(TT,'SetGlobalV','Channel',channel);
TT.SetEpocTimeFilterV('Freq',0,.2);
trs =  TT.GetValidTimeRangesV;
freq = invoke(TT, 'GetEpocsV', 'Freq', 0, 0, 10000);
levl = invoke(TT, 'GetEpocsV', 'Levl', 0, 0, 10000);
nlvl = invoke(TT, 'GetEpocsV', 'nlvl', 0, 0, 10000);


TT.SetGlobalV('RespectOffsetEpoc',0);
StartTime=-.25;     
StopTime=0;       
SpikesByTime=zeros(7,length(trs));
SpikesByTime(1,:)=levl(1,:);

TT.SetEpocTimeFilterV('Freq',-.25,.25);
trs = TT.GetValidTimeRangesV;
Chan1Snip = TT.ReadWavesOnTimeRangeV('Snip',1);    
for i=1:length(Chan1Snip(1,:))
    SpikesByTime(3,i)=length(find(Chan1Snip(:,i)));

end

SpikesByTime(5,:)=SpikesByTime(3,:);
SpikesByTime(7,:)=SpikesByTime(3,:);
SpikesByTime(2,:)=abs(StartTime);
            %Now that i have the baseline calculated I want to 
            %calculate the mean starting point for the new line
            

StartTime1=.0050;
StepSize1=.005;
StopTime1=StartTime1+StepSize1*4;
X1val=(StartTime1+StopTime1)/2;
SpikesByTime(4,:)=X1val+SpikesByTime(2,:);
k=0;
Summation=zeros(1,length(SpikesByTime(1,:)));
for j=StartTime1+StepSize1:StepSize1:StopTime1;   
    TT.SetEpocTimeFilterV('Freq',0,j);
    trs = TT.GetValidTimeRangesV;
    Chan1Snip = TT.ReadWavesOnTimeRangeV('Snip',1);
    for i=1:length(Chan1Snip(1,:))
        Summation(i)=Summation(i)+length(find(Chan1Snip(:,i)));
    end
    k=k+1;
end
Summation=Summation/(k+1);
SpikesByTime(5,:)=SpikesByTime(5,:)+Summation;  %This finds the mean Y value at


StartTime2=.01;
StepSize2=.005;
StopTime2=StartTime2+StepSize2*4;
X2val=(StartTime2+StopTime2)/2;
SpikesByTime(6,:)=X2val+SpikesByTime(2,:);
k=0;
Summation=zeros(1,length(SpikesByTime(1,:)));
for j=(StartTime2+StepSize2):StepSize2:StopTime2   
    TT.SetEpocTimeFilterV('Freq',0,j);
    trs = TT.GetValidTimeRangesV;
    Chan1Snip = TT.ReadWavesOnTimeRangeV('Snip',1);
    for i=1:length(Chan1Snip(1,:))
        Summation(i)=Summation(i)+length(find(Chan1Snip(:,i)));
    end
    k=k+1;
end
Summation=Summation/(k+1);
SpikesByTime(7,:)=SpikesByTime(7,:)+Summation;
close(figurea)


[Sorted,ix]=sort(SpikesByTime(1,:));
Sorted=SpikesByTime(:,ix);


    
OnsetDelay=zeros(3,length(Sorted));
OnsetDelay(1,:)=Sorted(1,:);
xvals=(0:.00001:.450);
y1vals=zeros(length(Sorted),length(xvals));
y2vals=zeros(length(Sorted),length(xvals));
for i=1:length(Sorted)       
    y1vals(i,:)=FindYvals(0,Sorted(2,i),0,Sorted(3,i),xvals);
    y2vals(i,:)=FindYvals(Sorted(4,i),Sorted(6,i),Sorted(5,i),Sorted(7,i),xvals);
    hold on
    [OnsetValue,Yval]=FindIntercept(xvals,y1vals(i,:),y2vals(i,:));
    OnsetDelay(2,i)=OnsetValue;
    OnsetDelay(3,i)=Yval;
end
OnsetDelay(2,:)=OnsetDelay(2,:)-.25;

for i=1:length(OnsetDelay(2,:))
    if OnsetDelay(2,i)<0 || OnsetDelay(2,i)>.05
        OnsetDelay(2,i)=NaN;
    end
end

end


function [yvals]=FindYvals(x1,x2,y1,y2,xvals)

    m=(y2-y1)/(x2-x1);
    b=y2-m*x2;
    yvals=m*xvals+b;
end

function [OnsetDelay,Yval]=FindIntercept(xvals,y1vals,y2vals)
y=y2vals-y1vals;
[Intercept,Index]=min(abs(y));
OnsetDelay=xvals(Index);
Yval=y2vals(Index);
end

