% Choice probability function 
% INPUTS: arrays of go/no-go resps (corr and incorr.) and
% their corresponding tone levels from Ram lab tone detection task (e.g. Dylla
% et al. 2013)
% OUTPUTS: choice probs, grand cps (via z-score), significance value
%
%
function [ChoiceProbabilities , Grand_ChoiceProbability, p, grandcpboot] = CP(a,b,tlvls)
 %hold off
[new_tlvlsC,new_tlvlsI] = Reshaped(a,b); %Finds data where there are >= 3 trials for each of Correct and Incorrect and reshapes data for use in ROC/CP

cList = new_tlvlsC(1,:); %List of tone levels where there are >= 3 trials for each of Correct and Incorrect

Correct = new_tlvlsC(2:end,:); %Takes away the first row which is the corresponding tone level, and leaves the spike data for Correct and Incorrect Trials
Incorrect = new_tlvlsI(2:end,:);

respPref = Correct'; 
respNonPref = Incorrect';


% Find the Max number of trials for Both Correct and Incorrect Responses 
counterC = zeros(1,length(new_tlvlsC(1,:)));
for j = 1:length(new_tlvlsC(1,:))
    for i = 2:length(new_tlvlsC(:,j))
        if isnan(new_tlvlsC(i,j)) == 0
            counterC(j) = counterC(j) + 1; % num corr
        end    
    end
end
counterI = zeros(1,length(new_tlvlsI(1,:)));
for j = 1:length(new_tlvlsI(1,:))
    for i = 2:length(new_tlvlsI(:,j))
        if isnan(new_tlvlsI(i,j)) == 0
            counterI(j) = counterI(j) + 1; % num incorr
        end    
    end
end


%% ROC curves

if max(respPref)>max(respNonPref)
    critList = 0:max(respPref(:));   %list of criterion values
else
    critList = 0:max(respNonPref(:)); % in the case where respNonPref is higher than Pref
end

pHit = zeros(length(critList),length(cList));
pFA =  zeros(length(critList),length(cList));
boots = 10;
pHitboot = zeros(length(critList),length(cList));
pFAboot =  zeros(length(critList),length(cList));


for cNum = 1:length(cList)
    % bootstrap at each tone level (cList), skipping NaNs
    corrnanidx=find(isnan(respPref(cNum,:)));
    inconanidx=find(isnan(respNonPref(cNum,:)));
    CorrBoot(cNum,:) = bootstrp(boots,@mean,respPref(cNum,1:corrnanidx-1));
    IncoBoot(cNum,:) = bootstrp(boots,@mean,respNonPref(cNum,1:inconanidx-1));
    
    % ROC curve at each tone level (bootstr and regular)
    for critNum = 1:length(critList)
        nTrialsC = counterC(cNum);
        nTrialsI = counterI(cNum);
        criterion = critList(critNum);
        pHit(critNum,cNum) = sum(respPref(cNum,:)>criterion)/nTrialsC;
        pFA(critNum,cNum) = sum(respNonPref(cNum,:)>criterion)/nTrialsI;
        pHitboot(critNum,cNum) = sum(CorrBoot(cNum,:)>criterion)/boots;
        pFAboot(critNum,cNum) = sum(IncoBoot(cNum,:)>criterion)/boots;
    end
    
end



%% This is the CP list (AUC of ROC curves) for each tone level 
ChoiceProbabilities = zeros(3,length(cList)); 
ChoiceProbabilities(1,:) = cList;
for i=1:length(cList)
    ChoiceProbabilities(2,i) = -trapz(pFA(:,i),pHit(:,i)); % trapz integrates the AUC
    ChoiceProbabilities(3,i) = -trapz(pFAboot(:,i),pHitboot(:,i));
end

%% Now lets do the grand CP calculation for all tone lvls 
% Z score the data for each of the tone lvls for both Correct (respPref) and Incorrect (respNonpref),
%after z scoring columnwise for each tone levl, we can put these z scores into a total Correct and total Incorrect  


% z_Correct = reshape(z_spikes_tlvls_C,[],1);
% z_Incorrect = reshape(z_spikes_tlvls_I,[],1);
% 
% z_Correct = rmmissing(z_Correct);
% z_Incorrect = rmmissing(z_Incorrect);  

% mean Value for the Zero Voltage tone level (No Tone)
zero_voltage_spikes = rmmissing(tlvls(2:end,1));
norm_zero_voltage_spikes = normalize(zero_voltage_spikes);
mean_quiet = mean(zero_voltage_spikes);
std_quiet = std(zero_voltage_spikes); 


%take the z score for each tone level, so we need columnwise
%z scores for these instead of entire
%array z score, 
%this also provides each of the mu and sigma for each tone level for both 
%correct and incorrect distrubutions 

% z score to incorr to the no tone response, z scor corr to responses on incorrect trials
for k = 1:length(cList) 

        
            I_mu(k)=nanmean(Incorrect(:,k));
            I_sigma(k)=nanstd(Incorrect(:,k));
            z_spikes_tlvls_I(:,k)=(Incorrect(:,k)-mean_quiet)./std_quiet;
        
            C_mu(k)=nanmean(Correct(:,k));
            C_sigma(k)=nanstd(Correct(:,k));
            z_spikes_tlvls_C(:,k)=(Correct(:,k)-nanmean(Incorrect(:,k)))./I_sigma(k);
       
      
end

k=1;
for k = 1:length(cList) 

        
            I_muboot(k)=nanmean(IncoBoot(k,:));
            I_sigmaboot(k)=nanstd(IncoBoot(k,:));
            z_spikes_tlvls_Iboot(:,k)=(IncoBoot(k,:)-mean_quiet)./std_quiet;
        
            C_muboot(k)=nanmean(CorrBoot(k,:));
            C_sigmaboot(k)=nanstd(CorrBoot(k,:));
            z_spikes_tlvls_Cboot(:,k)=(CorrBoot(k,:)-nanmean(IncoBoot(k,:)))./I_sigmaboot(k);
       
      
end

% Now that z score is taken columnwise we add to total arrays for correct and
% incorrect

zC_tot = rmmissing(reshape(z_spikes_tlvls_C,[],1));
zI_tot = rmmissing(reshape(z_spikes_tlvls_I,[],1)); 

zCorrboot = rmmissing(reshape(z_spikes_tlvls_Cboot,[],1));
zIncoBoot = rmmissing(reshape(z_spikes_tlvls_Iboot,[],1)); 

% Zscored Histogram
% figure(2)
% clf
% subplot(1,2,1);
% histogram(zC_tot,40)
% hold on 
% histogram(zI_tot,40) 


%% ROC analysis for grand CP, and randperm determined CP, to assess stat. significance
[Grand_ChoiceProbability] = CP_ROC(zC_tot,zI_tot);

%
parfor permct = 1:1:boots
    r = randi([1 500],10,1); % select a few random trials
    grandcpboot(1,permct) = CP_ROC(zCorrboot(r),zIncoBoot(r));
end


% determine significance of CP by calculating confidence intervals of the
% permuted data re chance (Greg DeAngelis method)
x=grandcpboot(1,:);
p=95;
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI = CIFcn(x,95);

if mean(grandcpboot)>0.5
    overlap = min(CI)-0.5;
elseif mean(grandcpboot)<0.5
    overlap = 0.5-max(CI);
end

if overlap>0
    p=0.001;
else 
    p=0.1;
end


% plot to check what boostrapped data look like
% close all;
% histogram(bC)
% hold on 
% histogram(bI)
% figure(2)
% histogram(zC_tot)
% hold on 
% histogram(zI_tot)

%[p, observeddifference, effectsize] = permutationTest(bC,bI,1000);

end