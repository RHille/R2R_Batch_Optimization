function r = HermiteIP(f,g,x)

weight = (1/sqrt(2*pi))*exp(-x^2/2);

h = collect(simplify(f*g*weight),x);

r = double(int(h,x,-inf,inf));