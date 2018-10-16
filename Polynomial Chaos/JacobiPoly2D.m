function r = JacobiPoly2D(x,y,order,a1,b1,alpha_p1,beta_p1,a2,b2,alpha_p2,beta_p2)

order1D = floor((1/2)*(1+sqrt(8*order-7)));

psi_x = JacobiPoly(x,order1D,a1,b1,alpha_p1,beta_p1);
psi_y = JacobiPoly(y,order1D,a2,b2,alpha_p2,beta_p2);

for i=1:order1D
    for j=1:order1D
        k = ((i+j-1)*(i+j))/2 - j + 1;
        psi_temp(k) = simplify(expand(psi_x(i)*psi_y(j)));
    end
end

psi = psi_temp(1:order+1);

r = psi;