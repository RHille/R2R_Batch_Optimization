function r = genincr(index, limit)

dimension = length(index);

for i=1:dimension
    if (index(i) < limit(i))
        index(i) = index(i) + 1;
        break;
    else
        index(i) = 0;
    end 
end

r = index;