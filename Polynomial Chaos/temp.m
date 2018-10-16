% this file is just a sandbox for playing with the functions.

% clear;
% 
% syms x y
% order = 9;
% a1 = -1; b1 = 1;
% a2 = -1; b2 = 1;
% alpha1 = 2; beta1 = 5;
% alpha2 = 2; beta2 = 5;
% 
% psi = pcepoly('hermite', x,order);
% 
% for i=1:order
%     for j=1:order
%         c(i,j) = HermiteIP(x, psi(i)*psi(j));
%     end
% end
% 
% c
% clear;
% 
% dimension = 3
% order = 10
% 
% 
% check = 0;
% order1D = 0;
% while(check < order)
%     order1D = order1D + 1;
%     check = check + nmultichoosek(dimension,order1D);
% end
% 
% order1D;
% 
% for i=1:order1D
%     x(i) = i;
% end
% 
% tempa = ones(1,dimension);
% tempa(1) = order1D
% 
% polytens = reshape(x, tempa);
% size(polytens)
% for i=2:dimension
%     i
%     tempa = ones(1,order1D);
%     tempa(i) = order1D;
%     polytens = polytens*reshape(x, tempa)
% end

% clear;
% dimension = 2;
% order = 6;
% for i=1:dimension
%     x(i) = sym(strcat('x',num2str(i)), 'real');
% end
% 
% check = 0;
% order1D = 0;
% while(check < order)
%     order1D = order1D + 1;
%     check = check + nmultichoosek(dimension,order1D);
% end
% 
% for i=1:dimension
%     mypolys(i,:) = HermitePoly(x(i),order1D);
% end
% 
% mypolys
% 
% limit = ones(1, dimension)*order1D;
% index = zeros(1, dimension);
% base = ones(1, dimension);
% 
% while(any(limit-index))
%     
%     
%     newindex = ndind(index)+1;
%     pcepolys(newindex) = mypolys(1,index(1)+1);
%     for i=2:dimension
%         pcepolys(newindex) = pcepolys(newindex)*mypolys(i,index(i)+1);
%     end
%     pcepolys(newindex) = simplify(expand(pcepolys(newindex)));
%     index = genincr(index,limit);
%     
%     limit-index
%     any(limit-index)
%     
% end
% 
% pcepolys
% 
% r = pcepolys(1:order);
% 
% r
clear;

dimension = 2;
order = 10;

for i=1:dimension
    x(i) = sym(strcat('x',num2str(i)), 'real');
end

alpha = [2 3];

a = [0 0];
b = [1 2];
alpha_p = [2 3];
beta_p = [3 4];

name = 'laguerre';

r = pcepolymd(name,x,order,alpha)

for i=1:order
    for j=1:order
        ips(i,j) = pcepolymdIP(name,r(i),r(j),x,alpha);
    end
end

ips





