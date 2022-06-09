
%% Plotting Grand CP values for durations values 
clc;close all;
dur = [0,50:12.5:200]';
last = 14; % num of durations used in the time window analysis
boot=0; % choose to see bootstrapped data instead



% I hate myself
% for meanct = 1:1:length(grandcpboot)
%     for colct = 1:1:length(grandcpboot(1,:))
%         deltaboot(meanct,colct) = mean(grandcpboot{meanct,colct});
%     end
% end



%% all durations

if boot==0
    alpha_CP_IC = alphacp(1:35,:);
    alpha_CP_CN = alphacp(36:end,:);
    bravo_CP_IC = bravocp(1:35,:);
    bravo_CP_CN = bravocp(36:end,:);
    charlie_CP_IC = charliecp(1:18,:);
    charlie_CP_CN = charliecp(19:end,:);
    delta_CP_IC = deltacp(1:20,:);
    delta_CP_CN = deltacp(21:end,:);
else
    alpha_CP_IC = alphaboot(1:35,:);
    alpha_CP_CN = alphaboot(36:end,:);
    bravo_CP_IC = bravoboot(1:35,:);
    bravo_CP_CN = bravoboot(36:end,:);
    charlie_CP_IC = charlieboot(1:18,:);
    charlie_CP_CN = charlieboot(19:end,:);
    delta_CP_IC = deltaboot(1:20,:);
    delta_CP_CN = deltaboot(21:end,:);
end

%% combine across monkeys
allIC=vertcat(alpha_CP_IC,bravo_CP_IC,charlie_CP_IC,delta_CP_IC);
moreICidx=find(allIC(:,last)>0.5); %indices of positive cps
lessICidx=find(allIC(:,last)<0.5); % indices of negative cps

allCN=vertcat(alpha_CP_CN,bravo_CP_CN,charlie_CP_CN,delta_CP_CN);
moreCNidx=find(allCN(:,last)>0.5); %indices of positive cps
lessCNidx=find(allCN(:,last)<0.5); % indices of negative cps

%% average across neurons whose cps go up vs. go down
moreICmean=mean(allIC(moreICidx,:));
moreICmean(2,:)=std(allIC(moreICidx,:))/sqrt(length(moreICidx));
moreICmean(3,:)=std(allIC(moreICidx,:));

lessICmean=mean(allIC(lessICidx,:));
lessICmean(2,:)=std(allIC(lessICidx,:))/sqrt(length(lessICidx));
lessICmean(3,:)=std(allIC(lessICidx,:));

moreCNmean=mean(allCN(moreCNidx,:));
moreCNmean(2,:)=std(allIC(moreCNidx,:))/sqrt(length(moreCNidx));

lessCNmean=mean(allCN(lessCNidx,:));
lessCNmean(2,:)=std(allCN(lessCNidx,:))/sqrt(length(lessCNidx));



%% time to significance

% p values separated by brain region
alpha_P_IC = alphap(1:35,:);
alpha_P_CN = alphap(36:end,:);
bravo_P_IC = bravop(1:35,:);
bravo_P_CN = bravop(36:end,:);
charlie_P_IC = charliep(1:18,:);
charlie_P_CN = charliep(19:end,:);
delta_P_IC = deltap(1:20,:);
delta_P_CN = deltap(21:end,:);

% put them together
allICp=vertcat(alpha_P_IC,bravo_P_IC,charlie_P_IC,delta_P_IC);
allCNp=vertcat(alpha_P_CN,bravo_P_CN,charlie_P_CN,delta_P_CN);

% find first time bin where p is significant
ICsig=[];
for ICsigrow = 1:1:length(allICp)
    
    sigtmp=find(allICp(ICsigrow,:)<0.1);
    if ~isempty(sigtmp)
    ICsig(ICsigrow,1)=dur(min(sigtmp),1);
    end

end

CNsig=[];
for CNsigrow = 1:1:length(allCNp)
    
    sigtmp=find(allCNp(CNsigrow,:)<0.1);
    if ~isempty(sigtmp)
    CNsig(CNsigrow,1)=dur(min(sigtmp),1);
    end

end

%% plotting
figure
% CPs @ 200 ms
subplot(3,2,1)
histogram(allIC(:,14),'BinWidth',0.05)
ylim([0,18]);
title('CP - IC')
subplot(3,2,2)
histogram(allCN(:,14),'BinWidth',0.05)
ylim([0,18]);
title('CP - CN')

% CP over time
subplot(3,2,3)

plot(dur,moreICmean(1,:))
hold on;
plot(dur,lessICmean(1,:))
title('CP vs Time - IC')
ylim([0.3,0.7]);
subplot(3,2,4)
plot(dur,moreCNmean(1,:))
hold on;
plot(dur,lessCNmean(1,:))
title('CP vs Time - CN')
ylim([0.3,0.7]);

% time to significance
edges = dur;
subplot(3,2,5)
h=histogram(ICsig,0:12.5:200,'Normalization','pdf')
title('Time to Significance (ms) - IC')
subplot(3,2,6)
h1=histogram(CNsig,0:12.5:200,'Normalization','pdf')
title('Time to Significance (ms) - CN')



[h,p]=kstest2(allIC(:,14),allCN(:,14))
          

