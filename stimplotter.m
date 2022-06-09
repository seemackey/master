figure
subplot(1,4,1)
cdfplot(thresh(nonstimidx))
hold on
cdfplot(thresh(stimidx))
title('Threshold')

subplot(1,4,2)
cdfplot(dr(nonstimidx))
hold on
cdfplot(dr(stimidx))
title('Dyn Range')

subplot(1,4,3)
cdfplot(rtint(nonstimidx))
hold on
cdfplot(rtint(stimidx))
title('RT int')

subplot(1,4,4)
cdfplot(rtslope(nonstimidx))
hold on
cdfplot(rtslope(stimidx))
title('RT slope')