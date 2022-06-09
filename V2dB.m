function dB=V2dB(V)

PlayListHandle=getappdata(0,'PlayListHandle');
v0=getappdata(PlayListHandle,'v0');
dB = 20.*log10((V)./v0);
end