function [modfreq_trim] = trimmer(modfreq,dur_trim)
    for trim = 1:1:length(modfreq)

        % don't trim the bin size for unmod noise
        if modfreq(1,trim) == 0
            modfreq(3,trim) = dur_trim;
            trim = trim+1;
        end

        % this array will have 1) mod freq, 2) # cycles, and 3) dur_trimmed
        modfreq(2,trim) = dur_trim/(1000/modfreq(1,trim)); % num cycles before trimming
        ncyc = modfreq(2,trim);
        if (ncyc-floor(ncyc))<0.5
        modfreq(3,trim) = floor(dur_trim-(((ncyc)-floor(ncyc))*(1000/modfreq(1,trim))));
        else
        modfreq(3,trim) = floor(dur_trim+((1-((ncyc)-floor(ncyc)))*(1000/modfreq(1,trim))));  
    end
    
    modfreq_trim=modfreq;
    
end