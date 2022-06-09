function [] = SaveThreshData(MyTank,MyBlock,NeuroThresh,BehavThresh,...
NeuroSlope,BehavSlope, WeightedDI, Loc, NoiseLevel, Monkey,ModNoise,...
RSstd,SpikeSlope,SpikeR, SpikedBSlope,SpikedBR,BFused,BehavThreshFit,...
NeuroThreshFit, RSBstd,BootstrapDI,BootstrapDIstd)




NoiseLevel=num2str(round(NoiseLevel));
if ModNoise==0
    Dirrectory= ['C:\Users\Kimaya\Desktop\ThresholdData\Thresholds_',Monkey,'_',Loc,'_', NoiseLevel, 'dB_Noise.txt'];
    AllThreshData='C:\Users\Kimaya\Desktop\ThresholdData\Thresholds.txt';
elseif ModNoise==1
    Dirrectory= ['C:\Users\Kimaya\Desktop\ThresholdData\Thresholds_',Monkey,'_',Loc,'_Modulated', NoiseLevel, 'dB_Noise.txt'];
    AllThreshData='C:\Users\Kimaya\Desktop\ThresholdData\ModulatedNoiseThresholds.txt';
end
    
fid=fopen(Dirrectory, 'a');
fclose all;  
fid=fopen(Dirrectory, 'r');
Everything=textscan(fid, '%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
fclose all;
%Identification Data includes date, monkey, block, record location, BF, and
%Noise level
tanks=Everything{1};
blocks=Everything{2};
locs=Everything{3};
ibfused=Everything{4};
nlvl=Everything{5};
%Neural Threshold Data
nthresh=Everything{6};
nthreshfit=Everything{7};
iRSstd=Everything{8};
nslope=Everything{9};
%Behavior Threshold Data (add thresh std in future...)
bthresh=Everything{10};
bthreshfit=Everything{11};
bthreshstd=Everything{12};
bslope=Everything{13};

di=Everything{14};
%Spike count vs Rxn time slope and R value
iSpikeSlope=Everything{15};
iSpikeR=Everything{16};
%dB level vs Rxn time slope and R value
iSpikedBSlope=Everything{17};
iSpikedBR=Everything{18};

iBootstrapDI= Everything{19};
iBootstrapDIstd = Everything{20};

MyTank=MyTank(end-10:end);

flag=0;
index=0;
for i=1:length(Everything{1})
    if strcmp(tanks{i},MyTank)==1
        
        if strcmp(blocks{i},MyBlock)==1
            
            flag=1;
            
            ibfused(i)=BFused;
            
            nthresh(i)=NeuroThresh;
            nthreshfit(i)=NeuroThreshFit;
            nslope(i)=NeuroSlope;
            iRSstd(i)=RSstd;
            
            bthresh(i)=BehavThresh;
            bslope(i)=BehavSlope;
            bthreshfit(i)=BehavThreshFit;
            bthreshstd(i)=0;
           
            di(i)=WeightedDI;
          
            iSpikeSlope(i)=SpikeSlope;
            iSpikeR(i)=SpikeR;
            iSpikedBSlope(i)=SpikedBSlope;
            iSpikedBR(i)=SpikedBR;
            
            iBootstrapDI(i)=BootstrapDI;
            iBootstrapDIstd(i)=BootstrapDIstd;
            
            index=i;
        end
    end
end

NoiseLevel=str2num(NoiseLevel);
j=length(tanks)+1;
if flag==0
    
    tanks(j)={MyTank};
    blocks(j)={MyBlock};
    locs(j)={Loc};
    
    ibfused(j)=BFused;
    nlvl(j)=NoiseLevel;
    
    nthresh(j)=NeuroThresh;
    nslope(j)=NeuroSlope;
    iRSstd(j)=RSstd;
    nthreshfit(j)=NeuroThreshFit;
    
    bthresh(j)=BehavThresh;
    bslope(j)=BehavSlope;
    bthreshfit(j)=BehavThreshFit;
    bthreshstd(j)=0;
    
    
    di(j)=WeightedDI;
    
    
    iSpikeSlope(j)=SpikeSlope;
    iSpikeR(j)=SpikeR;
    iSpikedBSlope(j)=SpikedBSlope;
    iSpikedBR(j)=SpikedBR;
    
    iBootstrapDI(j)=BootstrapDI;
    iBootstrapDIstd(j)=BootstrapDIstd;
    
end

fid=fopen(Dirrectory, 'w');
for j=1:length(tanks)    
    fprintf(fid, '\n');   
    fprintf(fid, '%s\t', (tanks{j}));
    fprintf(fid, '%s\t', (blocks{j}));
    fprintf(fid, '%s\t', (locs{j}));
    fprintf(fid, '%d\t', (ibfused(j)));
    fprintf(fid, '%d\t', (nlvl(j)));
    fprintf(fid, '%d\t', (nthresh(j)));
    fprintf(fid, '%d\t', (nthreshfit(j)));
    fprintf(fid, '%d\t', (iRSstd(j)));
    fprintf(fid, '%d\t', (nslope(j)));
    fprintf(fid, '%d\t', (bthresh(j)));
    fprintf(fid, '%d\t', (bthreshfit(j)));
    fprintf(fid, '%d\t', (bthreshfit(j)));
    fprintf(fid, '%d\t', (bslope(j)));
    fprintf(fid, '%d\t', (di(j))); 
    fprintf(fid, '%d\t', (iSpikeSlope(j)));
    fprintf(fid, '%d\t', (iSpikeR(j)));
    fprintf(fid, '%d\t', (iSpikedBSlope(j)));
    fprintf(fid, '%d\t', (iSpikedBR(j)));    
    fprintf(fid, '%d\t', (iBootstrapDI(j)));
    fprintf(fid, '%d\r', (iBootstrapDIstd(j)));
end
fclose(fid);


fid=fopen(AllThreshData, 'a');
fclose all;  
fid=fopen(AllThreshData, 'r');
Everything=textscan(fid, '%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
fclose all;

tanks=Everything{1};
blocks=Everything{2};
locs=Everything{3};
ibfused=Everything{4};
nlvl=Everything{5};
%Neural Threshold Data
nthresh=Everything{6};
nthreshfit=Everything{7};
iRSstd=Everything{8};
nslope=Everything{9};
%Behavior Threshold Data (add thresh std in future...)
bthresh=Everything{10};
bthreshfit=Everything{11};
bthreshstd=Everything{12};
bslope=Everything{13};

di=Everything{14};
%Spike count vs Rxn time slope and R value
iSpikeSlope=Everything{15};
iSpikeR=Everything{16};
%dB level vs Rxn time slope and R value
iSpikedBSlope=Everything{17};
iSpikedBR=Everything{18};

iBootstrapDI= Everything{19};
iBootstrapDIstd = Everything{20};

MyTank=MyTank(end-10:end);

flag=0;
index=0;
for i=1:length(Everything{1})
    if strcmp(tanks{i},MyTank)==1        
        if strcmp(blocks{i},MyBlock)==1            
            flag=1;
            
            ibfused(i)=BFused;
            
            nthresh(i)=NeuroThresh;
            nthreshfit(i)=NeuroThreshFit;
            nslope(i)=NeuroSlope;
            iRSstd(i)=RSstd;
            
            bthresh(i)=BehavThresh;
            bslope(i)=BehavSlope;
            bthreshfit(i)=BehavThreshFit;
            bthreshstd(i)=0;
           
            di(i)=WeightedDI;
          
            iSpikeSlope(i)=SpikeSlope;
            iSpikeR(i)=SpikeR;
            iSpikedBSlope(i)=SpikedBSlope;
            iSpikedBR(i)=SpikedBR;
            
            iBootstrapDI(i)=BootstrapDI;
            iBootstrapDIstd(i)=BootstrapDIstd;
            
            index=i;
        end
    end
end

j=length(tanks)+1;
if flag==0
    
    tanks(j)={MyTank};
    blocks(j)={MyBlock};
    locs(j)={Loc};
    
    ibfused(j)=BFused;
    nlvl(j)=NoiseLevel;
    
    nthresh(j)=NeuroThresh;
    nslope(j)=NeuroSlope;
    iRSstd(j)=RSstd;
    nthreshfit(j)=NeuroThreshFit;
    
    bthresh(j)=BehavThresh;
    bslope(j)=BehavSlope;
    bthreshfit(j)=BehavThreshFit;
    bthreshstd(i)=0;
    
    di(j)=WeightedDI;
 
    iSpikeSlope(j)=SpikeSlope;
    iSpikeR(j)=SpikeR;
    iSpikedBSlope(j)=SpikedBSlope;
    iSpikedBR(j)=SpikedBR;
    
    iBootstrapDI(j)=BootstrapDI;
    iBootstrapDIstd(j)=BootstrapDIstd;
end

fid=fopen(AllThreshData, 'w');
for j=1:length(tanks)    
    fprintf(fid, '\n');     
    fprintf(fid, '%s\t', (tanks{j}));
    fprintf(fid, '%s\t', (blocks{j}));
    fprintf(fid, '%s\t', (locs{j}));
    fprintf(fid, '%d\t', (ibfused(j)));
    fprintf(fid, '%d\t', (nlvl(j)));
    fprintf(fid, '%d\t', (nthresh(j)));
    fprintf(fid, '%d\t', (nthreshfit(j)));
    fprintf(fid, '%d\t', (iRSstd(j)));
    fprintf(fid, '%d\t', (nslope(j)));
    fprintf(fid, '%d\t', (bthresh(j)));
    fprintf(fid, '%d\t', (bthreshfit(j)));
    fprintf(fid, '%d\t', (bthreshfit(j)));
    fprintf(fid, '%d\t', (bslope(j)));
    fprintf(fid, '%d\t', (di(j))); 
    fprintf(fid, '%d\t', (iSpikeSlope(j)));
    fprintf(fid, '%d\t', (iSpikeR(j)));
    fprintf(fid, '%d\t', (iSpikedBSlope(j)));
    fprintf(fid, '%d\t', (iSpikedBR(j)));
    fprintf(fid, '%d\t', (iBootstrapDI(j)));
    fprintf(fid, '%d\r', (iBootstrapDIstd(j)));
    
end
fclose(fid);





