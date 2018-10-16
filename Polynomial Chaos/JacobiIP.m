function r = JacobiIP(f,g,x,a,b,alpha_p,beta_p)

weight = ((b-x)^alpha_p*(x-a)^beta_p)/...
    ((b-a)^(alpha_p+beta_p+1)*beta(alpha_p+1,beta_p+1));

h = collect(simplify(f*g*weight),x);

r = double(int(h,x,a,b));