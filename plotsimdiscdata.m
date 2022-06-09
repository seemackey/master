function [hr,mf,rt] = plotsimdiscdata(simdata500ms)
%% plotting simulated discrimination data

%% separate by monkey
    %simdata500msarr=table2array(simdata500ms);
    dari=find(simdata500ms.subj_idx=='dario');
    isil=find(simdata500ms.subj_idx=='isildur');

    
    %% remove nogos from RTs
    nogoidx=simdata500ms.rt>0.999;
    simdata500ms.rt(nogoidx)=NaN;
    
%% hit rate calculation and organize rts by condition/subj

    for subct = 1:1:2
        if subct == 1
                sub = isil;
        else
                sub = dari;
        end
        for condct = 1:1:length(unique(simdata500ms.NMF2(sub)))
            
            mfs = unique(simdata500ms.NMF2(sub)); % different mfs
            mfcurr = mfs(condct,1); % the one we're on in this iteration
            mf(condct,subct) = mfcurr; % storing later for plotting
            
            responses_currsubj = simdata500ms.response(sub); % all resps of this subj
            responses_currsubj_idx = (find(simdata500ms.NMF2(sub)==mfcurr)); % indices of resps at curr mf

            responses = responses_currsubj(responses_currsubj_idx); % responses of this subj at curr mf
            
            rt_currsubj = simdata500ms.rt(sub); % all resps of this subj
            rt_currsubj_idx = (find(simdata500ms.NMF2(sub)==mfcurr)); % indices of resps at curr mf

            rt{condct,subct} = rt_currsubj(rt_currsubj_idx); % responses of this subj at curr mf
            
            hr(condct,subct) = length(find(responses==1))/length(responses); % calculate hit rate
      
        end
    end
    
    
    
%% plotting
    % pretty sure this is dario
    figure
    subplot(2,2,1)
    scatter(mf(:,1),hr(:,1))
    title('Isildur Hit Rate')
    subplot(2,2,2)
    cdfplot(rt{4,1})
    hold on
    cdfplot(rt{5,1})
    hold on
    cdfplot(rt{7,1})
    legend('28 Hz','36 Hz','84 Hz')
    title('Isildur RTs')
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
    % pretty sure this is isildur
    subplot(2,2,3)
    scatter(mf(:,2),hr(:,2))
    title('Dario Hit Rate')
    subplot(2,2,4)
    cdfplot(rt{4,2})
    hold on
    cdfplot(rt{5,2})
    hold on
    cdfplot(rt{7,2})
    legend('28 Hz','36 Hz','84 Hz')
    title('Dario RTs')
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
end