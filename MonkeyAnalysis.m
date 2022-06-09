function varargout = MonkeyAnalysis(varargin)
% MONKEYANALYSIS M-file for MonkeyAnalysis.fig
%      MONKEYANALYSIS, by itself, creates a new MONKEYANALYSIS or raises the existing
%      singleton*.
%
%      H = MONKEYANALYSIS returns the handle to a new MONKEYANALYSIS or the handle to
%      the existing singleton*.
%
%      MONKEYANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MONKEYANALYSIS.M with the given input arguments.
%
%      MONKEYANALYSIS('Property','Value',...) creates a new MONKEYANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MonkeyAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MonkeyAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MonkeyAnalysis

% Last Modified by GUIDE v2.5 02-Nov-2018 17:08:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MonkeyAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @MonkeyAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes on button press in PlotTogetherPlot.
function PlotTogetherPlot_Callback(hObject, eventdata, handles)
UseOld=get(handles.OldDateCheck, 'Value');
OldDate=get(handles.OldDate,'String');
MonkeySelect=get(handles.MonkeySelect, 'Value');
[MonkeyDataTank, MonkeyName] = GetMonkeyNames(MonkeySelect,UseOld,OldDate);

RecordHandle=get(handles.LocationSelect, 'Value');
RecordSite= GetRecordSite(RecordHandle);
TrialNum=str2num(get(handles.TrialNum, 'String'));
loopcount=length(TrialNum);
SaveThresh=get(handles.SaveThresh,'Value');
%Monkeyplottogether(MonkeyName,MonkeyDataTank,TrialNum,RecordSite)


for i=1:loopcount
    MonkeyNandLevel(MonkeyName,MonkeyDataTank,TrialNum(i),RecordSite,SaveThresh);
end
    % --- Executes on button press in ResponseBlockPlot.
end

function ResponseBlockPlot_Callback(hObject, eventdata, handles) %#ok<DEFNU>
Type=get(handles.Behavior,'Value');
UseOld=get(handles.OldDateCheck, 'Value');
OldDate=get(handles.OldDate,'String');
MonkeySelect=get(handles.MonkeySelect, 'Value');
[MonkeyDataTank, MonkeyName] = GetMonkeyNames(MonkeySelect,UseOld,OldDate);

RecordHandle=get(handles.LocationSelect, 'Value');
RecordSite= GetRecordSite(RecordHandle);
TrialNum=str2num(get(handles.TrialNum, 'String'));
%a= strcat( 'G:\', MonkeyName, ' data\NeuroBehavior\',...
%    RecordSite, '\ex', MonkeyDataTank, '\tank', MonkeyDataTank);

MonkeyRespBlockPlotter(MonkeyName,MonkeyDataTank,TrialNum,RecordSite,Type);
end
function ThreshPlotGo_Callback(hObject, eventdata, handles)
MonkeySelect=get(handles.MonkeySelect2, 'Value');
LowestNoise=get(handles.LowestNoiseSelect, 'Value');
RecordHandle=get(handles.LocationSelect2, 'Value');
Location= GetRecordSite(RecordHandle);
PlotThreshCorr(LowestNoise,Location, MonkeySelect)
end
 

% --- Executes on button press in ResponseMapPlot.
function ResponseMapPlot_Callback(hObject, eventdata, handles)

UseOld=get(handles.OldDateCheck, 'Value');
OldDate=get(handles.OldDate,'String');
MonkeySelect=get(handles.MonkeySelect, 'Value');
[MonkeyDataTank, MonkeyName] = GetMonkeyNames(MonkeySelect,UseOld,OldDate);
FitLineSelect=get(handles.FitLineSelect,'Value');
RecordHandle=get(handles.LocationSelect, 'Value');
RecordSite= GetRecordSite(RecordHandle);
TrialNum=str2num(get(handles.TrialNum, 'String'));
Date=MonkeyDataTank(1:6);
StringArray=['b' 'r' 'k' 'm' 'b'];
StringArray2=['+' '*' '.' '+' '*'];
LoopCount=length(TrialNum);

scrsz = get(0,'ScreenSize'); %ScreenSize is a four-element vector: [left, bottom, width, height]:
switch LoopCount
    case 1
        figure('Position',[1 scrsz(4)/2 scrsz(3)/2 (scrsz(4)/2-75)])
    case 2
        figure('Position',[1 scrsz(4)/3 scrsz(3)/2 (scrsz(4)/1.5)-75]) 
    otherwise
        figure('Position',[1 1 scrsz(3)/2 (scrsz(4)-75)])    
end

hold on

for i=1:LoopCount
    [Data1 Data2 StimLevel warnflag, baseline]= MonkeyRMapPlotter(MonkeyName,MonkeyDataTank,TrialNum(i),RecordSite);    
    Freq(i,:)=Data1;
    Count(i,:)=Data2;   
    if warnflag==0
        FitX=(Freq(i,:));
        FitY=(Count(i,:));
        realXind = find(~isnan(FitX));
        realX = FitX(realXind);
        realY = FitY(realXind);
        realYind= find(~isnan(realY));
        realX = realX(realYind);
        realY = realY(realYind);
        realX=rot90(realX);
        realY=rot90(realY);
        flag=0;
        switch FitLineSelect
            case 1; flag=1;
            case 2; fitline=fit(realX,realY,'gauss4'); clc;                
            case 3; fitline=fit(realX,realY,'poly9');       %fit line...
            case 4; fitline=fit(realX,realY,'fourier6');       %fit line...
        end
        SaveResponseMapData(Data1,Data2,baseline,MonkeyName,RecordSite,TrialNum(i),MonkeyDataTank)
        StringArray3=strcat(StringArray(i),StringArray2(i));
        h(i)= subplot(LoopCount,1,i);
        hold on
        plot(Freq(i,:),Count(i,:),StringArray3, 'MarkerSize',6);
        baseline=zeros(1,length(Freq(i,:)))+baseline;
        plot(sort(Freq(i,:)),baseline,'k-')
        if flag==0
            hold on
            plot(fitline,'k--')
        end
        hold off
        p = get(h(i), 'Position');
        p(3) = p(3) + 0.05;
        p(4) = p(4) +.07;
        p(2) = p(2) -.03;
        set(h(i), 'pos', p);
        set(gca,'xscale','log');
        set(gca,'xlim',[min(Freq(i,:)) max(Freq(i,:))]);
        StimLevel=num2str(StimLevel);
        legend(StimLevel,'Location','NorthOutside')
    end
end
hold off
end


% --- Executes on button press in TankPlotterPlot.
function TankPlotterPlot_Callback(hObject, eventdata, handles)
Type=get(handles.Behavior,'Value');
UseOld=get(handles.OldDateCheck, 'Value');
OldDate=get(handles.OldDate,'String');
MonkeySelect=get(handles.MonkeySelect, 'Value');
[MonkeyDataTank, MonkeyName] = GetMonkeyNames(MonkeySelect,UseOld,OldDate);

RecordHandle=get(handles.LocationSelect, 'Value');
RecordSite= GetRecordSite(RecordHandle);
TrialNum=str2num(get(handles.TrialNum, 'String'));
%a= strcat( 'G:\', MonkeyName, ' data\NeuroBehavior\',...
%   RecordSite, '\ex', MonkeyDataTank, '\tank', MonkeyDataTank);

MonkeyTankPlotter(MonkeyName,MonkeyDataTank,TrialNum,RecordSite,Type)
figure;
createFitWB
end


% --- Executes on button press in CloseAnalysis.
function CloseAnalysis_Callback(hObject, eventdata, handles)

close MonkeyAnalysis
end



function [DataTank, MonkeyName]= GetMonkeyNames(MonkeySelect,UseOld,OldDate)
switch MonkeySelect
    case 1
        MonkeyName='Alpha';
        FileExt='a';
    case 2
        MonkeyName='Bravo';
        FileExt='b';
    case 3
        MonkeyName='Charlie';
        FileExt='c';
    case 4
        MonkeyName='Delta';
        FileExt='d';
    case 5
        MonkeyName='Echo';
        FileExt='e';
    case 6
        MonkeyName='Gatsby';
        FileExt='g';
    case 7
        MonkeyName='Lima';
        FileExt='l';
    case 8
        MonkeyName='Bilbo';
        FileExt='b';
    case 9
        MonkeyName='Haldir';
        FileExt='h';
    otherwise
        
end

if UseOld==0
    DateString=GetDateString();
    DataTank=strcat(DateString, FileExt);
else
    OldDate= num2str(OldDate);
    DataTank=strcat(OldDate, FileExt);
end
end


function [out]= GetDateString()
DateArray=clock;
year=num2str(DateArray(1)-2000);
month=num2str(DateArray(2));
if length(month)==1
    month= strcat('0',month);
end
day=num2str(DateArray(3));
if length(day)==1
    day= strcat('0',day);
end
out=strcat(year,month,day);
end



function [RecordSite]= GetRecordSite(input)
switch input
    case 1
        RecordSite='CN';
    case 2
        RecordSite='IC';
end
end


% --- Executes on button press in CloseFigures.
function CloseFigures_Callback(hObject, eventdata, handles)
global SaveData    
cab
SaveData.Flag=0;    
    
end



% --- Executes on button press in CumSpikeGo.
function CumSpikeGo_Callback(hObject, eventdata, handles)

UseOld=get(handles.OldDateCheck, 'Value');
OldDate=get(handles.OldDate,'String');
MonkeySelect=get(handles.MonkeySelect, 'Value');
[MonkeyDataTank, MonkeyName] = GetMonkeyNames(MonkeySelect,UseOld,OldDate);

RecordHandle=get(handles.LocationSelect, 'Value');
RecordSite= GetRecordSite(RecordHandle);
TrialNum=str2num(get(handles.TrialNum, 'String'));
MonkeyRespByTime(MonkeyName,MonkeyDataTank,TrialNum,RecordSite);
end



% --- Executes just before MonkeyAnalysis is made visible.
function MonkeyAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
end

function varargout = MonkeyAnalysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
end

function figure1_DeleteFcn(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

end
