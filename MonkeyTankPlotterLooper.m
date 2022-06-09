%Tankplotterloop
clc;close all;clear all;
tic
%clear the text file we'll write to
fopen('DRnRTdata.txt','w'); 

%paths to all the blocks we want to analyze
paths={
     
    
    '/Volumes/ramlab2018/Luthien data/Behavior/ex211214/tank211214l/OurData-3';
    '/Volumes/ramlab2018/Luthien data/Behavior/ex211214/tank211214l/OurData-2';
    '/Volumes/ramlab2018/Luthien data/Behavior/ex211214/tank211214l/OurData-1';
    '/Volumes/ramlab2018/Luthien data/Behavior/ex211201/tank211201l/OurData-4';
    '/Volumes/ramlab2018/Luthien data/Behavior/ex211201/tank211201l/OurData-3';
    '/Volumes/ramlab2018/Luthien data/Behavior/ex211201/tank211201l/OurData-2';
    '/Volumes/ramlab2018/Luthien data/Behavior/ex211208/tank211208l/OurData-1';
    
    };


  
%loop time, super efficient. you are smart for doing this.
for PlotLoopCount = 1:1:size(paths,1)
    All = MonkeyTankPlotterOfflineFxn(paths{PlotLoopCount,:});
end
toc