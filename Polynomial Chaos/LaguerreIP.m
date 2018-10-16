function r = LaguerreIP(f,g,x,alpha)

weight = (x^alpha*exp(-x))/(gamma(alpha+1));

h = collect(simplify(f*g*weight),x);

r = double(int(h,x,0,inf));