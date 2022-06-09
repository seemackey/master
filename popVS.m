function [done] = popVS(binnedspc,sortedspikes,stimspc,int)
%% inputs
% binnedspc - spike counts binned across the whole recording
% sortedspikes - spike times
% stimspc - spike count on each trial ( redundant with binnespc)
%
%% outputs
% done - vector strength at each modulation frequency


%% spike count at each mf
avestimspc = zeros(29,1);

n = 0; 
onefreq = 2.^[1:(1/3):10];
onefreq(1, 29) = 0; 
onefreq = roundn(onefreq,-7); 
modfreq=[64	25.39841683	20.1587368	812.7493386	1024	2	203.1873347	101.5936673	256	12.69920842	40.3174736	6.349604208	3.174802104	2.5198421	5.0396842	10.0793684	16	322.5397888	4	32	50.79683366	80.63494719	8	0	128	406.3746693	161.2698944	512	645.0795775	161.2698944	6.349604208	32	0	512	64	5.0396842	40.3174736	1024	203.1873347	20.1587368	80.63494719	8	812.7493386	2	50.79683366	16	406.3746693	10.0793684	256	128	101.5936673	645.0795775	4	25.39841683	12.69920842	322.5397888	3.174802104	2.5198421	50.79683366	80.63494719	12.69920842	40.3174736	32	161.2698944	10.0793684	512	645.0795775	16	322.5397888	812.7493386	20.1587368	1024	5.0396842	6.349604208	128	256	8	0	4	2.5198421	25.39841683	406.3746693	64	3.174802104	2	101.5936673	203.1873347	50.79683366	128	5.0396842	256	64	10.0793684	645.0795775	20.1587368	80.63494719	512	40.3174736	12.69920842	16	2.5198421	6.349604208	2	101.5936673	8	406.3746693	0	25.39841683	161.2698944	322.5397888	203.1873347	32	4	812.7493386	1024	3.174802104];
modfreq=roundn(modfreq,-7);
freq=modfreq;
whichspike = zeros(29,5); 
whichspike(:,5) = onefreq';

levelspikes = zeros(29,4); % this will end up having # spikes at each mf
for i = 1:29
    col = 0; 
    for n = 1:116
        if modfreq(1, n) == onefreq(1,i)
            col = col + 1; 
            levelspikes(i,col) = stimspc(n, 1);
            whichspike(i,col) = n;
        end
    end
end
%   
avestimspc = nanmean(levelspikes'); 
avestimspcstd = std(levelspikes'); 

%% Vector strength

%First get timestamps lined up
spikesortbysec = [];
col = 1;
for i = 1:int:length(binnedspc)
    if i > 1
        sofar = sum(binnedspc(1, 1:i-1));
        spikesortbysec(1:binnedspc(1, i), col) = sortedspikes(sofar+1:sofar+binnedspc(1,i), 1) ;
        col = col + 1;
    else 
        spikesortbysec(1:binnedspc(1, i), col) = sortedspikes(1:binnedspc(1,i), 1);
        col = col + 1;
    end
end      

% then run this to find the VS at each frequency
onefreq = 2.^[1:(1/3):10];
onefreq(1, 29) = 0; 
onefreq = roundn(onefreq,-7); 
modfreq=[64	25.39841683	20.1587368	812.7493386	1024	2	203.1873347	101.5936673	256	12.69920842	40.3174736	6.349604208	3.174802104	2.5198421	5.0396842	10.0793684	16	322.5397888	4	32	50.79683366	80.63494719	8	0	128	406.3746693	161.2698944	512	645.0795775	161.2698944	6.349604208	32	0	512	64	5.0396842	40.3174736	1024	203.1873347	20.1587368	80.63494719	8	812.7493386	2	50.79683366	16	406.3746693	10.0793684	256	128	101.5936673	645.0795775	4	25.39841683	12.69920842	322.5397888	3.174802104	2.5198421	50.79683366	80.63494719	12.69920842	40.3174736	32	161.2698944	10.0793684	512	645.0795775	16	322.5397888	812.7493386	20.1587368	1024	5.0396842	6.349604208	128	256	8	0	4	2.5198421	25.39841683	406.3746693	64	3.174802104	2	101.5936673	203.1873347	50.79683366	128	5.0396842	256	64	10.0793684	645.0795775	20.1587368	80.63494719	512	40.3174736	12.69920842	16	2.5198421	6.349604208	2	101.5936673	8	406.3746693	0	25.39841683	161.2698944	322.5397888	203.1873347	32	4	812.7493386	1024	3.174802104];
modfreq=roundn(modfreq,-7);
freq=modfreq;


r = zeros(116,1); 
radpersec = zeros(1, length(freq));


for Z = 1:length(freq)
    radpersec(1, Z) = (2*pi)*freq(1,Z);
    nZ=stimspc(Z,1);
    reps=spikesortbysec(1:nZ,Z);
    phaseangleZ=radpersec(1, Z) .* reps;
    xZ=cos(phaseangleZ);
    yZ=sin(phaseangleZ);
    r(Z,1)=sqrt((sum(xZ))^2+(sum(yZ))^2)/nZ;
end

avvs=zeros(29,1);
done= [];
for i = 1:29
    for n = 1:116
        if onefreq(1,i) == modfreq(1, n)
            avvs(i,1) = avvs(i,1) + r(n, 1); 
        end
    end
    done(i,1) = avvs(i,1)/4; % need to make this so the 4 isn't hard coded
end
end