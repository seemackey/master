% loop through the samsmodlong pooling fxn to get neurometric sensitivity
% of the whole population


%% file paths
% close all
% clear all
% clc
% 
% paths = {
%     
%      
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150611/tank150611c/OurData-7';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150612/tank150612c/OurData-15';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-2';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-3';
%    '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-11';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150618/tank150618c/OurData-17';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150724/tank150724c/OurData-8';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex150807/tank150807c/OurData-5';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151125/tank151125c/OurData-7';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151125/tank151125c/OurData-11';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151218/tank151218c/OurData-13';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151106/tank151106c/OurData-18';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex180523/tank180523c/OurData-5';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181005/tank181005c/OurData-22';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181009/tank181009c/OurData-13';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181012/tank181012c/OurData-16';
%    '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181016/tank181016c/OurData-23';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181018/tank181018c/OurData-31';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181113/tank181113c/OurData-14';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-9';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-18';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-22';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181129/tank181129c/OurData-5';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181130/tank181130c/OurData-6';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-13';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181120/tank181120c/OurData-17';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-5';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-10';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-19';
%     '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181205/tank181205c/OurData-23';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex181207/tank181207c/OurData-4';
%      '/Volumes/ramlab/Charlie data/Neurobehavior/CN/ex151218/tank151218c/OurData-13';
%     
%      };
%     
% 
% 
% %     
% % % import
% for imp = 1:1:size(paths)
%     data{imp,1} = TDTbin2mat(paths{imp,1});
% end


%% optional - load data
% load ICalldata.mat
% load ICpaths.mat
% load CNalldata.mat
% load CNpaths.mat
%% permutation loop
perms=1000;
thresh=zeros(perms,perms);
ct=0;
i=1;
VS = 0; % if we're doing population Vector Strength, not working yet

parfor i = 1:1:perms
    
    [pc] = samsmodlong_pooler_fxn(paths,data,VS);
    
    if max(pc(:,1)) > 0.759 % threshold from increase in resp
        ct = ct+1;
        thresh(i,1)=ThreshCalc(pc');
        thresh(i,1)=abs(thresh(i,1)-pc(1,2));
        
    elseif min(pc(:,1)) < 0.241 % threshold from decrease in resp
        ct = ct+1;
        pc(:,1)=1-pc(:,1);
        thresh(i,1)=ThreshCalc(pc');
        thresh(i,1)=abs(thresh(i,1)-pc(1,2));
    
    end
end

% proportion of neurons with a threshold
propthresh=ct/i

% find the ones with thresholds
nonzeroidx=find(thresh(:,1)>0);
meanthresh=mean(thresh(nonzeroidx,1))
medianthresh = median(thresh(nonzeroidx,1))

%calculate confidence intervals of the mean, using percentile method
x=thresh(nonzeroidx,1);
p=95;
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI = CIFcn(x,95)

% check distribution of data
histogram(thresh(:,1),'BinWidth',2)
beep on; beep;