clear all; 
%% First you have to import all the data
location = 'IC';
date = '150507';
MonkeyName = 'Delta';
monkey = 'd';
block = '8';

myTank = ['Z:\' MonkeyName ' data\Neurobehavior\' location '\ex' date '\tank' date monkey];
myBlock = ['OurData-' block];

% myTank = ['Z:\Delta data\Neurobehavior\IC\ex150518\tank150518d'];
% myBlock = ['OurData-16'];

% myTank = 'Z:\Echo data\Neurobehavior\IC\ex150519\tank150519e';
% myBlock = 'OurData-11';


% This depends on which version of TDT2mat you're using. Format it right. 
data = TDT2mat(myTank,myBlock,'T1',0,'T2',40,'TYPE',3); 

%data = TDT2mat(myTank,myBlock,0,40);

%% Now you want to see how many spikes there were over the whole window. 

allspikes = data.snips.eNeu.ts;

bins = 0.5:40.5; %This centers the bins over the right one second interval
figure; 
hist(allspikes,bins); %Let's look at it, then save that list
binnedspc = hist(allspikes,bins);

stimspc = zeros(20,1); %This is just the bins that have the stimulus in it. 
n = 0;
for i = 1:2:40;
    n = n+1;
    stimspc(n,1) = binnedspc(1,i);
end

avestimspc = zeros(10,1);

n = 0; 

for i = 1:10
    n = n + 1; 
    start = 2*i-1;
    stop = 2*i;
    avestimspc(n,1) = sum(stimspc(start:stop, 1))./2; 
end
    
    
modfreq=[32;32;256;256;512;512;4;4;2;2;64;64;16;16;128;128;8;8;0;0;];
freq=modfreq';

onefreq = [32 256 512 4 2 64 16 128 8 0];
y = avestimspc(10, 1); 

for i = 1:10 
    SSave(1, 1:i) = y;
end

figure; %Just looking at the bins that have stimulus played. 
semilogx(onefreq,avestimspc, '*r', onefreq, SSave, '-');


title('Spike Count');
xlabel('Modulation Frequency'); 
ylabel('# of Spikes per second')

%% Vector Strength! 

%Vector strength timing analysis script.
%for the current application, this will take the modulated noise
%signal and calculate the vector strength of the spikes from physiology
%according to the methods in Goldberg et al 1969.
%theta=phase angle
%x=cos(theta)
%y=sin(theta)
%direction is a measure of mean phase relation between the stimulus and
%discharge given by: phi=arctan(sum(y)/sum(x))+k*(pi), where k is zero or
%one depending on signs of sum(x) and sum(y)
%vector strength, r, is the length of the mean vector.
%r=squareroot(sum(x)^2+sum(y)^2)/n
%n=number of vectors


%First get timestamps lined up
spikesortbysec = [];
col = 1;
for i = 1:2:40
    if i > 1
        sofar = sum(binnedspc(1, 1:i-1));
        spikesortbysec(1:binnedspc(1, i), col) = allspikes(sofar+1:sofar+binnedspc(1,i), 1) ;
        col = col + 1;
    else 
        spikesortbysec(1:binnedspc(1, i), col) = allspikes(1:binnedspc(1,i), 1);
        col = col + 1;
    end
end

% collapse frequency
spikebyfreq = []; 
n = 1;
for i = 1:2:19 
    spikebyfreq(1:sum(stimspc(i:(i+1), 1)), n) = vertcat(spikesortbysec(1:stimspc(i, 1), i), spikesortbysec(1:stimspc((i+1), 1), (i+1)));
    n = n + 1; 
end
        

% then run this to find the VS at each frequency
r = zeros(10,1); 
radpersec = zeros(1, length(onefreq));

for Z = 1:length(onefreq)
    radpersec(1, Z) = (2*pi)*onefreq(1,Z);
    nZ=sum(stimspc(Z:(Z+1),1));
    reps=spikebyfreq(1:nZ,Z);
    phaseangleZ=radpersec(1, Z) .* reps;
    xZ=cos(phaseangleZ);
    yZ=sin(phaseangleZ);
    r(Z,1)=sqrt((sum(xZ))^2+(sum(yZ))^2)/nZ;
end

figure; 
semilogx(onefreq, r, 's');
title('Vector Strength'); 
xlabel('Modulation Frequency');
ylabel('Vector Strength'); 

%% Get Data for Spreadsheet. 
% This puts them in order of modfreq 
 AASpikeCount = zeros(1, 10);
 AAVS = zeros(1,9); 
AASpikeCount(1,1) = avestimspc(5, 1);
AASpikeCount(1,2) = avestimspc(4, 1);
AASpikeCount(1,3) = avestimspc(9, 1);
AASpikeCount(1,4) = avestimspc(7, 1);
AASpikeCount(1,5) = avestimspc(1, 1);
AASpikeCount(1,6) = avestimspc(6, 1);
AASpikeCount(1,7) = avestimspc(8, 1);
AASpikeCount(1,8) = avestimspc(2, 1);
AASpikeCount(1,9) = avestimspc(3, 1);
AASpikeCount(1,10) = avestimspc(10, 1);

AAVS(1,1) = r(5, 1);
AAVS(1,2) = r(4, 1);
AAVS(1,3) = r(9, 1);
AAVS(1,4) = r(7, 1);
AAVS(1,5) = r(1, 1);
AAVS(1,6) = r(6, 1);
AAVS(1,7) = r(8, 1);
AAVS(1,8) = r(2, 1);
AAVS(1,9) = r(3, 1);

