%% cleaning discrimination data for lba model

function All_clean = CleanDiscData(All,coh)
    % set up a structure, All is unclean, All_clean is clean
    
    if coh==0
    All_clean={};
    All_clean.cond=[];
    All_clean.correct=[];
    All_clean.rt=[];
    All_clean.stim=[];
    All_clean.response=[];
    
    % rank conditions (mod freqs) and assign to cond and stim
    [uniqueValues, ia, ic]=unique(All(:,3));
    
    All_clean.stim=ic;
    All_clean.cond(1:1:length(All),1)=ic;
    
    % code resps as corr or incorr
    for ct = 1:1:length(All)
        if All(ct,4)==2 && All(ct,3)==20 % corr rej
           All_clean.correct(ct,1)=1
        elseif All(ct,4)==2 && All(ct,3)>20 % Miss
            All_clean.correct(ct,1)=0
        elseif All(ct,4)==1 && All(ct,3)==20 % FA
            All_clean.correct(ct,1)=0
        elseif All(ct,4)==1 && All(ct,3)>20
            All_clean.correct(ct,1)=1 
        end
        
    end
    
    % rts
    All_clean.rt=(All(:,5)-All(:,6))*1000;
    All_clean.rt=round(All_clean.rt,0);
    
    %resp
    All_clean.response=All(:,4);
    
    
    %if we're not doing this in %coh do this part
    else
        
    All_clean={};
    All_clean.cond=[];
    All_clean.correct=[];
    All_clean.rt=[];
    All_clean.stim=[];
    All_clean.response=[];
    
    % rank conditions (mod freqs) and assign to cond and stim
    [uniqueValues, ia, ic]=unique(All(:,3));
    
    All_clean.stim=ic;
    All_clean.cond(1:1:length(All),1)=1;
    
    subd=All(:,3)-All(1,1);
    All_clean.coherence=(subd/max(subd)*100);
    
    % code resps as corr or incorr
    for ct = 1:1:length(All)
        if All(ct,4)==2 && All(ct,3)==20 % corr rej
           All_clean.correct(ct,1)=1
        elseif All(ct,4)==2 && All(ct,3)>20 % Miss
            All_clean.correct(ct,1)=0
        elseif All(ct,4)==1 && All(ct,3)==20 % FA
            All_clean.correct(ct,1)=0
        elseif All(ct,4)==1 && All(ct,3)>20
            All_clean.correct(ct,1)=1 
        end
        
    end
    
    % rts
    All_clean.rt=(All(:,5)-All(:,6))*1000;
    All_clean.rt=round(All_clean.rt,0);
    
    %resp
    All_clean.response=All(:,4);    
  
    end
    
    % want to get rid of catch trials?
%      catchidx=find(All_clean.cond==1);
%      All_clean.rt(catchidx,:)=[];
%      All_clean.response(catchidx,:)=[];
%      All_clean.stim(catchidx,:)=[];
%      All_clean.correct(catchidx,:)=[];
%      All_clean.cond(catchidx,:)=[];

end
