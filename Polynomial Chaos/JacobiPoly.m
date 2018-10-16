function r = JacobiPoly(x,order,a,b,alpha_p,beta_p)

psi(1) = sym(1);
for i=2:order
    psi(i) = simplify(expand(...
        ((-1)^(i-1))/(2^(i-1)*factorial(i-1))...
        * diff((b-x)^((i-1)+alpha_p)*(x-a)^((i-1)+beta_p),i-1)...
        * (b-x)^(-alpha_p)*(x-a)^(-beta_p)));
end

r = psi;