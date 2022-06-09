
function [stimspc,sortedspikes] = mfSort(data,dur,int)
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



%% changing length of analysis window based on mod. freq
modfreq=[64	25.39841683	20.1587368	812.7493386	1024	2	203.1873347	101.5936673	256	12.69920842	40.3174736	6.349604208	3.174802104	2.5198421	5.0396842	10.0793684	16	322.5397888	4	32	50.79683366	80.63494719	8	0	128	406.3746693	161.2698944	512	645.0795775	161.2698944	6.349604208	32	0	512	64	5.0396842	40.3174736	1024	203.1873347	20.1587368	80.63494719	8	812.7493386	2	50.79683366	16	406.3746693	10.0793684	256	128	101.5936673	645.0795775	4	25.39841683	12.69920842	322.5397888	3.174802104	2.5198421	50.79683366	80.63494719	12.69920842	40.3174736	32	161.2698944	10.0793684	512	645.0795775	16	322.5397888	812.7493386	20.1587368	1024	5.0396842	6.349604208	128	256	8	0	4	2.5198421	25.39841683	406.3746693	64	3.174802104	2	101.5936673	203.1873347	50.79683366	128	5.0396842	256	64	10.0793684	645.0795775	20.1587368	80.63494719	512	40.3174736	12.69920842	16	2.5198421	6.349604208	2	101.5936673	8	406.3746693	0	25.39841683	161.2698944	322.5397888	203.1873347	32	4	812.7493386	1024	3.174802104];
modfreq=roundn(modfreq,-7);

%calculates time-windows (modfreq(trim(3,:)) so no resps to incomp. AM are used
% [modfreq_trim] = trimmer(modfreq,dur*1000);
% 
% modfreq_trim(3,:)=modfreq_trim(3,:)/1000;

stimspc=[];
% counting spikes
% for ts_ct = 1:1:length(modfreq_trim)
%    stimspc(ts_ct,1) = length(find(sortedspikes(:,1)<=((ts_ct-1)+modfreq_trim(3,ts_ct))&(sortedspikes(:,1))>ts_ct-1));
% end

for ts_ct = 1:1:length(modfreq)
   stimspc(ts_ct,1) = length(find(sortedspikes(:,1)<=((ts_ct-1)+0.5)&(sortedspikes(:,1))>ts_ct-1));
end
end