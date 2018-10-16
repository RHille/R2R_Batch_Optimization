function r = hermitepolymd(x, order)

dimension = length(x);

check = 0;
order1D = 0;
while(check < order)
    order1D = order1D + 1;
    check = check + nmultichoosek(dimension,order1D);
end

for i=1:dimension
    mypolys(i,:) = HermitePoly(x(i),order1D);
end

limit = ones(1, dimension)*order1D;
index = zeros(1, dimension);
base = ones(1, dimension);

while(any(limit-index))
    newindex = ndind(index)+1;
    pcepolys(newindex) = mypolys(1,index(1)+1);
    for i=2:dimension
        pcepolys(newindex) = pcepolys(newindex)*mypolys(i,index(i)+1);
    end
    simplify(expand(pcepolys(newindex)));
    index = genincr(index,limit);
    
end

r = pcepolys(1:order);


