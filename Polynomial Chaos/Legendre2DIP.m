function r = Legendre2DIP(f,g,x,y)

h = expand(simplify(f*g));

r = double(int(int(h,x,-1,1),y,-1,1));