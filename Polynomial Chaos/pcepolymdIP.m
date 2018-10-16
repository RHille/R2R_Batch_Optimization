function r = pcepolymdIP(name, varargin)

f = varargin{1};
g = varargin{2};
x = varargin{3};

dimension = length(x);

name = lower(name);
weight = sym('1');
switch(name)
    case 'legendre'
        llim = -1*ones(1,dimension); ulim = 1*ones(1,dimension);
        for i=1:dimension
            weight = weight*sym('1/2');
        end
    case 'hermite'
        llim = -inf*ones(1,dimension); ulim = inf*ones(1,dimension);
        for i=1:dimension
            weight = weight*(1/sqrt(2*pi))*exp(-x(i)^2/2);
        end
    case 'laguerre'
        llim = 0*ones(1,dimension); ulim = inf*ones(1,dimension);
        alpha = varargin{4};
        for i=1:dimension
            weight = weight*((x(i)^(vpa(alpha(i)))*exp(-x(i)))/(gamma(alpha(i)+1)));
        end
    case 'jacobi'
        a = varargin{4};
        b = varargin{5};
        alpha_p = varargin{6};
        beta_p = varargin{7};
        
        llim = a; ulim = b;
        
        for i=1:dimension
            weight = weight*((b(i)-x(i))^alpha_p(i)*(x(i)-a(i))^beta_p(i))/...
                ((b(i)-a(i))^(alpha_p(i)+beta_p(i)+1)*beta(alpha_p(i)+1,beta_p(i)+1));
        end
    otherwise
        error('Invalid argument.');
end

%weight = vpa(simplify(expand(weight)));
%weight = vpa(weight);

h = f*g*weight;

if(size(findsym(h)) ~= [0 0])
    h = collect(expand(h));
end

myint = int(h,x(1),llim(1),ulim(1));
for i=2:dimension
    myint = int(myint,x(i),llim(i),ulim(i));
end


r = double(myint);










