function r = LegendreIP(f,g,x)

h = collect(simplify(f*g),x);

r = double(int(h,x,-1,1));