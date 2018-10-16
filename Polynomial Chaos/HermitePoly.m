function r = HermitePoly(x,order)

psi(1) = sym(1);
for i=2:order
    psi(i) = simplify(expand(...
        ((-1)^(i-1))...
        * diff(exp(-x^2/2),i-1)...
        * exp(x^2/2)));
end

r = psi;