function r = pcepolymd(name, varargin)

x = varargin{1};
order = varargin{2};

dimension = length(x);

check = 0;
order1D = 0;
while(check < order)
    order1D = order1D + 1;
    check = check + nmultichoosek(dimension,order1D);
end

name = lower(name);
switch(name)
    case 'legendre'
        for i=1:dimension
            mypolys(i,:) = LegendrePoly(x(i),order1D);
        end
    case 'hermite'
        for i=1:dimension
            mypolys(i,:) = HermitePoly(x(i),order1D);
        end
    case 'laguerre'
        alpha = varargin{3};
        for i=1:dimension
            mypolys(i,:) = LaguerrePoly(x(i),order1D,alpha(i));
        end
    case 'jacobi'
        a = varargin{3};
        b = varargin{4};
        alpha_p = varargin{5};
        beta_p = varargin{6};
        for i=1:dimension
            mypolys(i,:) = JacobiPoly(x(i),order1D,a(i),b(i),alpha_p(i),beta_p(i));
        end
    otherwise
        error('Invalid argument.');
end


limit = ones(1, dimension)*(order1D-1);
index = zeros(1, dimension);

while(any(limit-index))
    newindex = ndind(index)+1;
    pcepolys(newindex) = mypolys(1,index(1)+1);
    for i=2:dimension
        pcepolys(newindex) = pcepolys(newindex)*mypolys(i,index(i)+1);
    end
    pcepolys(newindex) = simplify(expand(pcepolys(newindex)));
    index = genincr(index,limit);
    
end
for i=1:order
    pcepolys(i) = vpa(pcepolys(i));
end

r = pcepolys(1:order);


