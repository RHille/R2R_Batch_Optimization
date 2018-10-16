function r = LaguerrePoly(x,order,alpha)

psi(1) = sym(1);

for i=2:order
    psi(i) = simplify(expand(...
        (1/factorial(i-1))*diff(exp(-x)*x^((i-1)+alpha),i-1)*exp(x)*x^(-alpha)...
        ));
end

r = psi;