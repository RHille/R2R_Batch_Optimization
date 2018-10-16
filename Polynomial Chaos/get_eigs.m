function [omegas,eigvs] = get_eigs(m,a,c)

tol = 1e-10;
epsilon = 1e-6;


zips = zeros(1,m);

f_even = @(x) c - x*tan(a*x);
f_odd = @(x) x + c*tan(a*x);

for i=0:2:m-1
    xL = (pi/(2*a))*(i)+epsilon;
    xR = (pi/(2*a))*(i+1)-epsilon;
    
    zips(i+1) = bisect(xL,xR,f_even);
    
    xL = (pi/(2*a))*(i+1)+epsilon;
    xR = (pi/(2*a))*(i+2)-epsilon;
    
    zips(i+2) = bisect(xL,xR,f_odd);

end
omegas = zips;

eigv = zeros(1,m);
for i=1:m
    eigv(i) = (2*c)/(zips(i)^2 + c^2);
end
eigvs = eigv;