
function [binnedspc,stimspc,sortedspikes] = mfSort(data,dur,int)
%% SORT SPIKES

% To get just the ones from the newly sorted file. If you don't want the
% sorted one, you can change the sortcode number to 0 or you can change all
% of the places where it says sortedspikes back to allspikes below this. 

% old way
allspikes = data.snips.eNeu.ts;
% n=1;
% sortedspikes = []; 
% for i = 1:length(data{1,1}.snips.eNeu.sortcode)
%     if data{1,1}.snips.eNeu.sortcode(i,1) == 0 % 3 because it's sortcode is under the number 3
%         sortedspikes(n,1) = allspikes(i,1);
%         n = n+1;
%     end
% end

% new way -CM
sortidx=find(data.snips.eNeu.sortcode==0);
sortedspikes(sortidx,1)=allspikes(sortidx,1);



%% old way of binning %%
% bins = 0.25:0.5:116.25; %This centers the bins over the right one second interval
% figure; 
% hist(sortedspikes,bins); %Let's look at it, then save that list
% binnedspc = hist(sortedspikes,bins);
%
%% new way of binning where duration can be manipulated (uses bin edges not center)
% figure
edges=0:dur:116.5; % left bin edges
binnedspc=histogram(sortedspikes,edges);

% now try to manipulate bin width/duration of window where we count spikes
% figure
% histogram(sortedspikes,'BinWidth',0.5,'BinEdges',edges)

% assign spikes to this variable
binnedspc=binnedspc.Values;

%% COUNT SPIKES
stimspc = zeros(116,1); %This is just the bins that have the stimulus in it. 
n = 0;
for i = 1:int:length(binnedspc)
    n = n+1;
    stimspc(n,1) = binnedspc(1,i);
end
end