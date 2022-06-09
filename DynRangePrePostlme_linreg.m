%% 
datatable=table(dr,dur,exposed,freq,monkey);

lme3=fitlme(datatable,'dr~dur*exposed+freq*exposed+(1|monkey)')

% now the lme is made, but to show people what's going on we're going to
% make a linear regression, which can be plotted

linreg=fitlm(datatable,'dr~dur*exposed+freq*exposed')

% plot the residuals
figure
plotResiduals(lme,'fitted')

% get idx of drs that were pre vs post exp.
datatable_arr=table2array(datatable);
preidx=find(datatable_arr(:,3)==0);
postidx=find(datatable_arr(:,3)==1);

preDRs=datatable_arr(preidx,1); %allocate
postDRs=datatable_arr(postidx,1); %allocate

% plot the interactions with data on top
plotInteraction(linreg,'exposed','freq','predictions')
hold on
scatter(datatable_arr(preidx,4),preDRs)
hold on
scatter(datatable_arr(postidx,4),postDRs)

% plot the interactions with data on top
figure
plotInteraction(linreg,'exposed','dur','predictions')
hold on
scatter(datatable_arr(preidx,2),preDRs)
hold on
scatter(datatable_arr(postidx,2),postDRs)

% plot 3d plots so we can see pre/post of specific durations/frequencies
figure
plot3(dur(preidx),preDRs,freq(preidx),'o')
hold on
plot3(dur(postidx),postDRs,freq(postidx),'x')
legend('pre','post')
grid on

% plot pre/post dyn range at diff durations
figure
scatter(dur(preidx),preDRs)
hold on
scatter(dur(postidx),postDRs)

% plot pre/post dyn range at diff durations
figure
scatter(dur(preidx),mean(preDRs))
hold on
scatter(dur(postidx),mean(postDRs))