
%% Plotting Grand CP values for durations values 


A = load("D:\Jackson\Chase Stuff\Alpha_CP_Analysis\Dur = 12.5 -200, 12.5 Step, 0.0 window\Alpha_Grand_CPs.mat");%35 IC
B = load("D:\Jackson\Chase Stuff\Bravo_CP_Analysis\Dur = 12.5 -200, 12.5 Step, 0.0 window\Bravo_Grand_CPs.mat");%35 IC
C = load("D:\Jackson\Chase Stuff\Charlie_CP_Analysis\Dur = 12.5 -200, 12.5 Step, 0.0 window\Charlie_Grand_CPs.mat");%18 IC
D = load("D:\Jackson\Chase Stuff\Delta_CP_Anlaysis\Dur = 12.5 -200, 12.5 Step, 0.0 window\Delta_Grand_CPs.mat");%20 IC



% Need to pull out all IC vs CN tanks
alpha_CP_IC = A.Grand_CP(1:35,:);
alpha_CP_CN = A.Grand_CP(36:end,:);
bravo_CP_IC = B.Grand_CP(1:35,:);
bravo_CP_CN = B.Grand_CP(36:end,:);
charlie_CP_IC = C.Grand_CP(1:18,:);
charlie_CP_CN = C.Grand_CP(19:end,:);
delta_CP_IC = D.Grand_CP(1:20,:);
delta_CP_CN = D.Grand_CP(21:end,:);

x1 = 12.5:12.5:200; %Duration values to match the anaylsis files

a_y = mean(cell2mat(alpha_CP_IC));
b_y = mean(cell2mat(bravo_CP_IC));
c_y = mean(cell2mat(charlie_CP_IC));
d_y = mean(cell2mat(delta_CP_IC));

avg  = median([a_y;b_y;c_y;d_y]);
%     plot(x1,a_y,x1,b_y,x1,c_y,x1,d_y)
    plot(x1,avg)
%     legend('alpha','bravo','charlie','delta')
% for plot_ct = 1:n_ct %Loop through each of the valid neuron tanks
%    ;
%     plot(x1,y)
%     hold on
% end