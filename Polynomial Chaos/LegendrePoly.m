function r = LegendrePoly(x,order)

psi(1) = sym(1);
if(order > 1)
    psi(2) = x;
    for n=2:order-1
        psi(n+1) = simplify(expand(...
            (1/(n))*((2*(n-1)+1)*x*psi(n) - (n-1)*psi(n-1))...
            ));
    end
end

r = psi;