% Volt level to dBSPL
% Written by Andy Hrnicek under guidance of P.I. Ram Ramachandran, PhD.
% Fall 2010
% 
% This function does a conversion from volts to decible SPL for the tone
% levels, it could be done really simply to add the noise conversion in
% here if you wanted to do a signal to noise ratio (Tlvl in dBSPL - Nlvl in
% dBSPL).
% 
%%% INPUTS
% ToneLevels - array of tone level voltages used in trial in ascending
% order
% alpha - length of ToneLevels array
%%% OUTPUTS
% dbSPL - array of tone leves in ascending order in dBSPL 

function dbSPL = dBtoSPL(ToneLevels, alpha)
%ToneLevels;
dbSPL = zeros(1,alpha);
%v0 = .79/(10^(74/20)); % On Hylebos-PC in room 040
v0=.79/(10^(82/20)); % On Yakima in Room 034
for i = 1:alpha
    dbSPL(i) = 20*log10((ToneLevels(i))/v0);
     if dbSPL(i)==-Inf
         dbSPL(i)=-40;
     end
end
dbSPL=dbSPL;
 %dbSPL=dbSPL+46;
% disp('WRONG LVLS')
end