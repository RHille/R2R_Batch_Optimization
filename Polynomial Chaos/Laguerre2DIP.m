function r = Laguerre2DIP(f,g,x,y,alpha1,alpha2)

weight1 = (x^alpha1*exp(-x))/(gamma(alpha1+1));
weight2 = (y^alpha2*exp(-y))/(gamma(alpha2+1));

h = expand(simplify(f*g*weight1*weight2));

r = double(int(int(h,x,0,inf),y,0,inf));

