% pulls out number of nonzero values in an array

function arraylength = nonzerolength(array)

for i = length(array):-1:1
    if array(i) == 0
        array(i) = [];
    else
        break
    end
end

arraylength = length(array);

end