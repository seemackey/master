function [] =OutputImportantData(HitRateData,PCData,Responses,levl,freq,threshold,RTmin,RTmax)

HitRateData(1,:)=PCData(2,:);
HitRateVsLevel=rot90(HitRateData([1,5],2:end));
PCData(:,1)=[];
PCData=rot90(flipud(PCData));
HitRateVsLevel;
PCData;       



for i=1:length(levl(1,:))   
    TrueCases=find(Responses(2,:)>levl(2,i));
    if ~isempty(TrueCases)
        levl(3,i)=Responses(1,TrueCases(1));
        levl(4,i)=Responses(2,TrueCases(1))-levl(2,i);    
    end
end


levl(1,:)=round(dBtoSPL(levl(1,:),length(levl(1,:))));

[levl(1,:),Indeces]=sort(levl(1,:));
levl(2:end,:)=levl(2:end,Indeces);
a=find(levl(3,:)==1);
CorrectResponses=levl([1,4],a);
ReactionMode=mode(levl(4,:));

DeleteThese=find(CorrectResponses(2,:)>ReactionMode);
CorrectResponses(:,DeleteThese)=[];

deletethis=find(CorrectResponses(1,:)<(threshold-4));
CorrectResponses(:,deletethis)=[];
% 
CorrectResponses(2,:)=CorrectResponses(2,:)*1000;
  %figure

% PLOT RT VS LVL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % scatter(CorrectResponses(1,:),CorrectResponses(2,:))
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% allocate RTs in a range of lvls%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    RT1=[];
% % % 
% %    min=56; % collect RTs starting above this lvl
% %    max=58; %collect RTs below this lvl
% % % % % % % MEAN OR MEDIAN RT AT TWO LEVELs
% % % % % % %this loop allocates RTs by a given level so we can calc mean/median
%  for aa=1:length(CorrectResponses)
%     if CorrectResponses(1,aa)>RTmin&&CorrectResponses(1,aa)<RTmax %specify the lower and upper limit for the level
%         RT1(1,aa)=CorrectResponses(2,aa);
%     else
%         continue
%     end
% end
%  RT1(:,~any(RT1,1))=[];   %% deletes zeros
%  RT1=RT1'
%  
%  reps=length(RT1) % printing this out to make sure we're getting enough reps, probably don't use less than 12-15
%  
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%%%%%%%%%%%%%%%%%%%%% end of RT allocation %%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   
% %   
% %   
% %   
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% %   %%% plotting RT cdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  fontSize = 20;
% % 
% % % Compute the histogram of RT
%  [countsA, binsA] = hist(RT1);
% % 
% % % Compute the cumulative distribution function of RT
%  cdfA = cumsum(countsA) / sum(countsA);
% % 
% % % Plot the probability distribution of RTs
%  figure
%  subplot(2,2, 1);
%  bar(binsA, countsA);
%  title('Histogram of RT', 'FontSize', fontSize);
%  ylabel('Count RT', 'FontSize', fontSize);
%  xlabel('Values of RT', 'FontSize', fontSize);
%  grid on;
% % 
% % % Plot the cumulative distribution function of RTs
%  subplot(2,2, 3);
% % 
%  [h,stats]=cdfplot(RT1); 
% % 
% % % putting the cdf data in the workspace so it can be exported to figures
% h;
% yvals=get(h,'YData');
% xvals=get(h,'XData');
% xvals=xvals';
% yvals=yvals';
% assignin('base','yvals',yvals);
% assignin('base','xvals',xvals);
% % 
% % % fitting the RTs with a weibull function 
% createFitWB_RT(xvals,yvals,mean(RT1))
% stats
% title('CDF of RT', 'FontSize', fontSize);
% grid on;
% ylabel('Percentage RT (/100)', 'FontSize', fontSize);
% xlabel('Values of RT', 'FontSize', fontSize);
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%% end of RT cdf plotting %%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%   
%   
% % check to make sure enough reps
% reps= length(RT1)
% min
% max


fclose all


levl;

% RT VS LEVEL FIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=polyfit(CorrectResponses(1,:),CorrectResponses(2,:),1);

r = p(1) .* CorrectResponses(1,:) + p(2);
slope=p(1)
intercept=p(2)

% settings for the fit

    
% 
% hold on
% plot (CorrectResponses(1,:),r)
fid1 = fopen('DRnRTdata.txt','a+');
fprintf(fid1,'%4.4f\t',slope);
fclose(fid1);

fid1 = fopen('DRnRTdata.txt','a+');
fprintf(fid1,'%4.4f\t',intercept);
fclose(fid1);


end
