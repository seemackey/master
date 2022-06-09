% Written by Andy Hrnicek under guidance of P.I. Ram Ramachandran, PhD.
% Fall 2010
% 
% this code finds all of the different tone levels used in the block being
% looked at and sorts them into ascending order. it returns an array of
% these levels and the number of them in ToneLevel and alpha respectively.
%
%%% INPUTS
% levl - tone level by trial
%%% OUTPUTS
% ToneLevel - tone levels in ascending order in volts
% alpha - number of tone levels used in this block

function [ToneLevel alpha] = TLevlFindNOrder(levl)

if strcmp(num2str(levl(1,1)),'NaN')
    error('Datatank returned no values');
end

tonelevels=[]; %creates array
ToneLevel=sort(levl(1,:)); %this becomes an infinite loop without the sort...
while ~isempty(ToneLevel) %iterates through until ToneLevel is and empty array
tonelevels=[min(ToneLevel),tonelevels];
    while(ToneLevel(1)==tonelevels(1)) %deletes all instances of a value added to tonelevels
    ToneLevel(1)=[];
        if isempty(ToneLevel) %prevents error due to exceeding matrix dimensions
            break
        end
    end
end

for i=1:length(tonelevels) %reorders them for ease of use
    ToneLevel=[tonelevels(i) ToneLevel];
end

alpha = length(tonelevels);

end