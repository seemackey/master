function [ pos,neg ] = SplitData( data )

[A,B] = size(data);

posIter = 1;
negIter = 1;

for ii = 1:A
    if data(ii,1) > 0
        pos(posIter,1) = data(ii,1);
        pos(posIter,2) = data(ii,2);
        posIter = posIter+1;
    elseif data(ii,1) < 0
        neg(negIter,1) = data(ii,1);
        neg(negIter,2) = data(ii,2);
        negIter = negIter+1;
    end
end

end

