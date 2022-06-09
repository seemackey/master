function [ReSampledMean, ReSampledStd]= MonkeyBStrapResamp(BinnedArrayIn,tlvls)

[Height,Width]=size(BinnedArrayIn);
SpikeOccCount=zeros(Height-1,1);
TotalSamples=1000;
v0 = .79/(10^(74/20));




for i=2:Height
    SpikeOccCount(i-1)=sum(BinnedArrayIn(i,:));
end

BFArray=zeros(length(tlvls)-1,TotalSamples+1);
BFArray(:,1)=tlvls(1,2:end);
for j=1:TotalSamples
    CurrentData=zeros(Height,Width);
    CurrentData(1,:)=BinnedArrayIn(1,:);
    for i=2:Height
        for h=1:SpikeOccCount(i-1)
            RandValue=ceil(rand*SpikeOccCount(i-1));
            count=0;
            g=1;
            while (count<RandValue)
                count=count+BinnedArrayIn(i,g);
                if count>=RandValue
                    CurrentData(i,g)=CurrentData(i,g)+1;
                end
                g=g+1;
            end
        end
        
    end
    BFArray(:,j+1)=ROCAnal(CurrentData,tlvls);
end



MeanAndStd=zeros(length(BFArray(:,1)),3);
MeanAndStd(:,1)=BFArray(:,1);
for i=1:length(tlvls)-1    
    MeanAndStd(i,1)= 20*log10((MeanAndStd(i,1))/v0);
    BFArray(i,1)=20*log10((BFArray(i,1))/v0);
    MeanAndStd(i,2)=mean(BFArray(i,2:end));
    MeanAndStd(i,3)=std(BFArray(i,2:end));
end

Thresh=zeros(1,TotalSamples);
Slopes=Thresh;
for i=1:TotalSamples
    [XVal,~,Slope]=FindThreshold(BFArray(:,1),BFArray(:,i+1));
    Thresh(i)=XVal;
    Slopes(i)=Slope;
end

ReSampledMean=nanmean(Thresh);
ReSampledStd=nanstd(Thresh);






function output=ROCAnal(CurrentData,tlvls)

[Height,Width]=size(CurrentData);
ROCarray=zeros(Height-1,Width);

for j=1:Height-1
    for i=1:Width
        d=sum(CurrentData(j+1,1:i));
        e=sum(CurrentData(j+1,:));
        ROCarray(j,i)=1-(d/(e+eps));
    end
end

lengthROC=length(ROCarray(1,:));
heightROC=length(ROCarray(:,1));
sums=zeros(heightROC-1,2);
for j=2:heightROC
    for i=(lengthROC):-1:2
        X1=ROCarray(1,i);
        X2=ROCarray(1,i-1);
        Y2=ROCarray(j,i-1);
        Y1=ROCarray(j,i);
        deltaX=X2-X1;
        deltaY=Y2-Y1;
        sums(j-1,2)=sums(j-1,2)+Y1*deltaX+deltaX*deltaY/2;
        sums(j-1,1)=tlvls(1,j);
    end    
end
output=sums(:,2);


function [XVal,YVal, slope]=FindThreshold(XDataIn, YDataIn)
threshold=.76;
a = length(YDataIn);
YDifference=zeros(1,a-1);


for i=2:a
    if YDataIn(i-1)<.76 && YDataIn(i)>=.76
        YDifference(i-1)=YDataIn(i)-YDataIn(i-1);
    else
        YDifference(i-1)=NaN;
    end
end


[Value,Indices]=max(YDifference);
slope=Value/(XDataIn(Indices+1)-XDataIn(Indices));
b=YDataIn(Indices)-slope*XDataIn(Indices);
xvalues=(XDataIn(Indices):.1:XDataIn(Indices+1));
yvalues=slope*xvalues+b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           To find the closest value to threshold, what you can do is subtract .76
%                           from all numbers, find the absolute value of all numbers, and then find
%                           the minimum index value in the array. This will find the closest
%                           intercept of a non continuous function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,index]=min(abs(yvalues-threshold));
XVal=xvalues(index);
XVal=(XVal*slope)/slope;
YVal=yvalues(index);
YVal=(YVal*slope)/slope;

