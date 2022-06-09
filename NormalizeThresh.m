function [normthresh] = NormalizeThresh(thresholds)
normthresh=zeros(length(thresholds),2);
normthreshct=1;
    for normthreshct = 1:7:length(thresholds)
        normthresh(normthreshct:normthreshct+6,1)=thresholds(normthreshct:normthreshct+6,1)-thresholds(normthreshct,1);
        normthresh(normthreshct:normthreshct+6,2)=[200;100;50;25;12.5;6.50;3.25];
        if length(normthreshct)==108
            break
        end
    end
end