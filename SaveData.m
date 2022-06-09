function SaveData(SaveData)

Date=GetDateString();
MonkeyName=SaveData.MonkeyName;
Test=SaveData.TestType;
ToneLevel=num2str(SaveData.ToneLevel);
CF=num2str(SaveData.CF);
ToneModulation=num2str(SaveData.ToneMod);
NoiseType=SaveData.NoiseType;
NoiseLevel=num2str(SaveData.NoiseLevel);
NoiseMod=num2str(SaveData.ModNoiseFreq);
Type=SaveData.DataType;
Comments1=get(handles.Comments,'String');
FinalString=[MonkeyName,'	',Date,'	',Block,'	',Test,'	',ToneLevel,'	',CF,'	',ToneModulation,...
    '	',NoiseType,'	',NoiseLevel,'	',BehavOrNeuro,'	',NoiseMod,'	',Comments1];

Dirrectory= 'C:\Users\Kimaya\Desktop\Recording Info\RecordingInfo.txt';
fid=fopen(Dirrectory, 'a+');
fprintf(fid, '%s\r\n', FinalString);
fclose(fid);
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