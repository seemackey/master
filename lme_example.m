%% Variables%%%%
%%%%%% put your data into matlab as column variables %%%%%%%%%
% example:
% X=[3;4;5;6];
% Z=[1;1;1;1];

% put those variables into a table
Loc=str2num(char(Locs));
t = table(Monkey,cfs,Thresholds,N2threshold,slopes,N2slope,Loc,durs);

% run this line of code, you can name 'lme' whatever you want
% follow the formatting in the mathworks link for 'fitlme' 
lme = fitlme(t,'RTslope~asym+syn*freq+NW+syn*NW+(1|monkey)')