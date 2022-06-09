%% Loops through a function that runs TDT's basic data extraction function
clc;close all;clear all;
tic

%% initialize
All = zeros(1,1);
monkey=5; %alphabetical (Bi, Ch, De, Ga, Ha)
exp=1; %monkey was noise exposed 0,1

% paths to all the data
paths={'/Volumes/ramlab2018/Isildur data/Behavior/ex210616/tank210616i/OurData-10';
      '/Volumes/ramlab2018/Dario data/Behavior/ex210617/tank210617d/OurData-3';
      '/Volumes/ramlab2018/Dario data/Behavior/ex210611/tank210611d/OurData-2';
    }; 

%% loop that pulls out all trial data
check=0; % this tells the loop whether or not to concatenate arrays
for PlotLoopCount = 1:1:size(paths,1)
    All_tmp = RawExtract(paths{PlotLoopCount,:},monkey)';
    
    if check==0 
        All=All_tmp;
        check=check+1;
    else
    All=cat(1,All,All_tmp);
    end
end
All(:,3)=dBtoSPL(All(:,3),length(All(:,3))); %convert to dB, make sure it's set to the right room!

All(:,9)=exp; 
toc