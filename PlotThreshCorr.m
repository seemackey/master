           
function  [] = PlotThreshCorr(NoiseLevel,Location, Monkey)


fclose all;
switch Monkey
    case 1, Monkey='Alpha';        
    case 2, Monkey='Bravo';        
    case 3, Monkey='Both';
    case 4, Monkey='Charlie';    
    case 5, Monkey='Delta';
    case 6, Monkey='All';
end

ColorString='k*';

switch NoiseLevel
    case 1
        NoiseLevel='24';
        if strcmp(Monkey,'Both')==1
            [x y x2 y2]=SingleData(NoiseLevel,Location,'Alpha');
            [x3 y3 x4 y4]=SingleData(NoiseLevel,Location,'Bravo');
            x=vertcat(x,x3);
            y=vertcat(y,y3);
            x2=vertcat(x2,x4);
            y2=vertcat(y2,y4);
        else
            [x y x2 y2]=SingleData(NoiseLevel,Location,Monkey);
        end
    case 2
        if strcmp(Monkey,'Both')==1
            All=0;
            [x y x2 y2]=MultiData(Location,'Alpha',All);
            [x3 y3 x4 y4]=MultiData(Location,'Bravo',All);
            x=vertcat(x,x3);
            y=vertcat(y,y3);
            x2=vertcat(x2,x4);
            y2=vertcat(y2,y4);
        else
            All=1;
            [x y x2 y2]=MultiData(Location,Monkey,All);
        end
    case 3
        NoiseLevel='34';
        if strcmp(Monkey,'Both')==1
            [x y x2 y2]=SingleData(NoiseLevel,Location,'Alpha');
            [x3 y3 x4 y4]=SingleData(NoiseLevel,Location,'Bravo');
            x=vertcat(x,x3);
            y=vertcat(y,y3);
            x2=vertcat(x2,x4);
            y2=vertcat(y2,y4);
        else
            [x y x2 y2]=SingleData(NoiseLevel,Location,Monkey);
        end
    case 4
        NoiseLevel='44';
        if strcmp(Monkey,'Both')==1
            [x y x2 y2]=SingleData(NoiseLevel,Location,'Alpha');
            [x3 y3 x4 y4]=SingleData(NoiseLevel,Location,'Bravo');
            x=vertcat(x,x3);
            y=vertcat(y,y3);
            x2=vertcat(x2,x4);
            y2=vertcat(y2,y4);
        else
            [x y x2 y2]=SingleData(NoiseLevel,Location,Monkey);
        end
    case 5
        NoiseLevel='54';
        if strcmp(Monkey,'Both')==1
            [x y x2 y2]=SingleData(NoiseLevel,Location,'Alpha');
            [x3 y3 x4 y4]=SingleData(NoiseLevel,Location,'Bravo');
            x=vertcat(x,x3);
            y=vertcat(y,y3);
            x2=vertcat(x2,x4);
            y2=vertcat(y2,y4);
        else
            [x y x2 y2]=SingleData(NoiseLevel,Location,Monkey);
        end
    case 6
        if strcmp(Monkey,'Both')==1
            All=1;
            [x y x2 y2]=MultiData(Location,'Alpha',All);
            [x3 y3 x4 y4]=MultiData(Location,'Bravo',All);
            x=vertcat(x,x3);
            y=vertcat(y,y3);
            x2=vertcat(x2,x4);
            y2=vertcat(y2,y4);
        else
            All=1;
            [x y x2 y2]=MultiData(Location,Monkey, All);
            
        end

end



[XHisto,YHisto]=CorrelationHistogram(x,y);



%ThresholdPlot
figure
subplot(3,1,1)
hold on
plot(x,y,ColorString, x,x,'-r')
a1=polyfit(x,y,1);
fitline=a1(1)*x+a1(2);
plot(x,fitline,'-k')
a='Thresholds';
title(a ,... 
  'FontWeight','bold')
ylabel('Behavioral Thresholds');
xlabel('Neuro Thresholds');
hold off

subplot(3,1,2)
hold on


realXind = find(~isnan(x2));
realX = x2(realXind);
realY = y2(realXind);
realYind= find(~isnan(realY));
realX = realX(realYind);
realY = realY(realYind);
a=polyfit(realX,realY,1);

fitline=a(1)*x2+a(2);
corrcoef(x2,y2);
plot(x2,y2,ColorString)
plot(x2,fitline,'--')
a='Slopes';
title(a ,... 
  'FontWeight','bold')
ylabel('Behavioral Slopes');
xlabel('Neuro Slopes');
hold off

subplot(3,1,3)
hold on
bar(XHisto,YHisto)
hold off
ylabel('Occurences');
xlabel('Behavior Thresholds - Neural Thresholds');

function [x y x2 y2]=SingleData(NoiseLevel,Location,Monkey)

Dirrectory=['C:\Users\Kimaya\Desktop\ThresholdData\Thresholds_',...
    Monkey,'_',Location,'_', NoiseLevel, 'dB_Noise.txt'];
if exist(Dirrectory)==2
    fid=fopen(Dirrectory,'r');
    Everything=textscan(fid, '%s%s%s%f%f%f%f%f%f');
    fclose(fid);
else
    Everything={NaN NaN NaN NaN NaN NaN NaN NaN NaN};
end
x=Everything{5};
y=Everything{6};
x2=Everything{7};
y2=Everything{8};


function [x y x2 y2]= MultiData(Location,Monkey, All)

if All==1
    NoiseLevel='24';
    Dirrectory=['C:\Users\Kimaya\Desktop\ThresholdData\Thresholds_',...
        Monkey,'_',Location,'_', NoiseLevel, 'dB_Noise.txt'];
    if exist(Dirrectory)==2
        fid=fopen(Dirrectory,'r');
        Everything0=textscan(fid, '%s%s%s%f%f%f%f%f%f');
        fclose(fid);
    else
        Everything0={NaN NaN NaN NaN NaN NaN NaN NaN NaN};
    end
elseif All==0
    Everything0={NaN NaN NaN NaN NaN NaN NaN NaN NaN};
end

NoiseLevel='34';
Dirrectory=['C:\Users\Kimaya\Desktop\ThresholdData\Thresholds_',...
    Monkey,'_',Location,'_', NoiseLevel, 'dB_Noise.txt'];
if exist(Dirrectory)==2
    fid=fopen(Dirrectory,'r');
    Everything1=textscan(fid, '%s%s%s%f%f%f%f%f%f');
    fclose(fid);
else
    Everything1={NaN NaN NaN NaN NaN NaN NaN NaN NaN};
end

NoiseLevel='44';
Dirrectory=['C:\Users\Kimaya\Desktop\ThresholdData\Thresholds_',...
    Monkey,'_',Location,'_', NoiseLevel, 'dB_Noise.txt'];
if exist(Dirrectory)==2
    fid=fopen(Dirrectory,'r');
    Everything2=textscan(fid, '%s%s%s%f%f%f%f%f%f');
    fclose(fid);
else
    Everything2={NaN NaN NaN NaN NaN NaN NaN NaN NaN};
end

NoiseLevel='54';
Dirrectory=['C:\Users\Kimaya\Desktop\ThresholdData\Thresholds_',...
    Monkey,'_',Location,'_', NoiseLevel, 'dB_Noise.txt'];
if exist(Dirrectory)==2
    fid=fopen(Dirrectory,'r');
    Everything3=textscan(fid, '%s%s%s%f%f%f%f%f%f');
    fclose(fid);
else
    Everything3={NaN NaN NaN NaN NaN NaN NaN NaN NaN};
end



x=vertcat(Everything0{5},Everything1{5}, Everything2{5},Everything3{5});
y=vertcat(Everything0{6},Everything1{6}, Everything2{6},Everything3{6});
x2=vertcat(Everything0{7},Everything1{7}, Everything2{7},Everything3{7});
y2=vertcat(Everything0{8},Everything1{8}, Everything2{8},Everything3{8});

function [xval yval]= CorrelationHistogram(x,y)
b=.05; %Step Size
HistoX=(-1:b:1);
HistoY=zeros(1,length(HistoX));
PoolScale=max(max(x),max(y));
HistoX=PoolScale.*HistoX;
ScaledB=b*PoolScale;

for j=1:length(x)
    k=1;
    for i=HistoX(1):ScaledB:HistoX(end-1);
        if (y(j)-x(j)-i*ScaledB)>0 && (y(j)-x(j)-(i+ScaledB)*ScaledB)<=0
            HistoY(k)=HistoY(k)+1;
        end
        k=k+1;
    end
end
           
xval=HistoX;
adjustment=(xval(2)-xval(1))/2;
xval=xval+adjustment;
yval=HistoY;
xval;
            
            
