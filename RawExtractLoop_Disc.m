%% Loops through a function that extracts and sorts TDT data
clc;close all;clear all;
tic

%% initialize
All = zeros(1,1);
monkey=1; %alphabetical (Da, Is)


% paths to all the data
paths={'/Volumes/ramlab2018/Dario data/Behavior/ex210617/tank210617d/OurData-3';
       '/Volumes/ramlab2018/Dario data/Behavior/ex210618/tank210618d/OurData-2';
       %'/Volumes/ramlab2018/Dario data/Behavior/ex210506/tank210506d/OurData-3';
        %'/Volumes/ramlab2018/Dario data/Behavior/ex211114/tank211114d/OurData-3';
      '/Volumes/ramlab2018/Isildur data/Behavior/ex210616/tank210616i/OurData-10';
    }; 

%% loop that pulls out all trial data
check=0; % this tells the loop whether or not to concatenate arrays
for PlotLoopCount = 1:1:size(paths,1)
    All_tmp = RawExtract_Disc(paths{PlotLoopCount,:},monkey)';
    
    if check==0 
        All=All_tmp;
        check=check+1;
    else
    All=cat(1,All,All_tmp);
    end
end

toc