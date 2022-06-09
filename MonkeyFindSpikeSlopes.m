function [slope rval slope2 rval2]=MonkeyFindSpikeSlopes(Sorted,NeuroThresh)
 v0 = .79/(10^(82/20));



SortedPlotter=Sorted;
for i=1:length(Sorted)
    if Sorted(3,i)==2
        SortedPlotter(:,i)=NaN;
    end
    if Sorted(5,i)>.5
        SortedPlotter(:,i)=NaN;
    end
    if Sorted(5,i)<.2
        SortedPlotter(:,i)=NaN;
    end
end

for i=1:length(SortedPlotter)
    SortedPlotter(1,i) = 20*log10((SortedPlotter(1,i))/v0);
end

SortedPlotter=RemoveNans(SortedPlotter);
SortedPlotter=SortAscending(SortedPlotter);
ReSortedPlotter=SortedPlotter;
a=length(SortedPlotter);

for i=a:-1:1
    if SortedPlotter(1,i)<NeuroThresh(1)
        ReSortedPlotter(:,i)=[];
    end
end

p=polyfit(ReSortedPlotter(4,:),ReSortedPlotter(5,:),1);
r=p(1).*ReSortedPlotter(4,:)+p(2);
rval=corrcoef(ReSortedPlotter(4,:),ReSortedPlotter(5,:));
figure

p2=polyfit(ReSortedPlotter(1,:),ReSortedPlotter(5,:),1);
rval2=corrcoef(ReSortedPlotter(1,:),ReSortedPlotter(5,:));


hold on
plot(ReSortedPlotter(4,:),ReSortedPlotter(5,:),'k*')
plot(ReSortedPlotter(4,:),r,'-')
hold off

slope=p(1);
rval=rval(1,2);
slope2=p2(1);
rval2=rval2(1,2);
end

function x=RemoveNans(x) %and inf's from columns
x=rot90(x);
x(any(isnan(x),2),:) = [];
x(any(isinf(x),2),:) = [];
x=rot90(x,3);
end


function output=SortAscending(x)
[~,d2]=sort(round(x(1,:)));
output=x(:,d2);
end