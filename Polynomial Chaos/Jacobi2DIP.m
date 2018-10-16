function r = Jacobi2DIP(f,g,x,y,a1,b1,alpha_p1,beta_p1,a2,b2,alpha_p2,beta_p2)

weight1 = ((b1-x)^alpha_p1*(x-a1)^beta_p1)/...
    ((b1-a1)^(alpha_p1+beta_p1+1)*beta(alpha_p1+1,beta_p1+1));
weight2 = ((b2-y)^alpha_p2*(y-a2)^beta_p2)/...
    ((b2-a2)^(alpha_p2+beta_p2+1)*beta(alpha_p2+1,beta_p2+1));

h = expand(simplify(f*g*weight1*weight2));

r = double(int(int(h,x,a1,b1),y,a2,b2));