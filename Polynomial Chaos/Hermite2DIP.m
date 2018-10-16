function r = Hermite2DIP(f,g,x,y)

weight1 = (1/sqrt(2*pi))*exp(-x^2/2);
weight2 = (1/sqrt(2*pi))*exp(-y^2/2);

h = expand(simplify(f*g*weight1*weight2));

r = double(int(int(h,x,-inf,inf),y,-inf,inf));
