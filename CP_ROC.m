function [Grand_ChoiceProbability] = CP_ROC(zC_tot,zI_tot)

    all_z =  cat(1,zC_tot,zI_tot);  % Add all of the zscores into one array , C and I 
    
    z_critList = min(all_z):.0001:max(all_z);   %list of criterion values *** THIS NEEDS TO BE FIXED TO ACCOUNT FOR NEGATIVE VALUES
    
    z_pHit = zeros(length(z_critList),1);
    z_pFA =  zeros(length(z_critList),1);
    
        for critNum = 1:length(z_critList)
            nTrialsC = length(zC_tot);
            nTrialsI = length(zI_tot);
            criterion = z_critList(critNum);
            z_pHit(critNum) = sum((zC_tot(:))>criterion)/nTrialsC;
            z_pFA(critNum) = sum((zI_tot(:))>criterion)/nTrialsI;
        end
Grand_ChoiceProbability = -trapz(z_pFA,z_pHit);
end