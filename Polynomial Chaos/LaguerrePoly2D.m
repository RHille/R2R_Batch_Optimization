function r = LaguerrePoly2D(x,y,order,alpha1,alpha2)

order1D = floor((1/2)*(1+sqrt(8*order-7)));

psi_x = LaguerrePoly(x,alpha1,order1D);
psi_y = LaguerrePoly(y,alpha2,order1D);

for i=1:order1D
    for j=1:order1D
        k = ((i+j-1)*(i+j))/2 - j + 1;
        psi_temp(k) = simplify(expand(psi_x(i)*psi_y(j)));
    end
end

psi = psi_temp(1:order+1);

r = psi;