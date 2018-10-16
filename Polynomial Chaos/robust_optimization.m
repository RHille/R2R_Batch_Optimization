function [x,f,eflag,outpt,variance] = robust_optimization(U_history,new_corr,K,theta1_PCE,theta2_PCE,pci_func,pci_den,i_2,para_selection,Nom_values,sampling,bio_nom,sf_nom_mod)

% if nargin == 1 % No options supplied
%     opts = [];
% end

U0 = U_history(i_2,:);
if i_2 >= 2
    U0old = U_history(i_2-1,:);
else
    U0old = U0;
end

xLast = []; % Last place computeall was called
myf = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint

fun = @objfun; % the objective function, nested below
cfun = @constr; % the constraint function, nested below

gleq=[1	-0.978228658146	0.0556685671162
2	-0.887062599768	0.125580369465
3	-0.730152005574	0.186290210928
4	-0.519096129207	0.233193764592
5	-0.269543155952	0.26280454451
6	0	0.272925086778
7	0.269543155952	0.26280454451
8	0.519096129207	0.233193764592
9	0.730152005574	0.186290210928
10	0.887062599768	0.125580369465
11	0.978228658146	0.0556685671162];

epsi_1=gleq(:,2); epsi_2=gleq(:,2);
epsi1=[]; epsi2=[]; w=[];

for j=1:length(epsi_1)
    epsi1=[epsi1; ones(length(epsi_2),1)*epsi_1(j)];
    epsi2=[epsi2; gleq(1:end,2)];
    w=[w; gleq(j,3)*gleq(:,3)];
end
theta1=subs(theta1_PCE); theta2=subs(theta2_PCE);

% Call fmincon

lb=[0.001 0.001];
ub=[100 100];

options=optimset('Display','iter','Algorithm','active-set',...
    'TolX',1e-4,'TolFun',1e-4,'TolCon',1e-4,'AlwaysHonorConstraints','bounds');
[x,f,eflag,outpt] = fmincon(fun,U0,[],[],[],[],lb,ub,cfun,options,Nom_values,new_corr,K,theta1,theta2,epsi1,epsi2,w,pci_func,pci_den,para_selection,sampling,bio_nom,sf_nom_mod);

function y = objfun(x,Nom_values,new_corr,K,theta1,theta2,epsi1,epsi2,w,pci_func,pci_den,para_selection,sampling,bio_nom,sf_nom_mod)
    if ~isequal(x,xLast) % Check if computation is necessary
        [myf,myc,myceq,S,variance] = compute_robust_loop_input_uncertainty(x,Nom_values,new_corr,K,theta1,theta2,epsi1,epsi2,w,pci_func,pci_den,para_selection,sampling,bio_nom,sf_nom_mod);
        xLast = x;
    end
    % Now compute objective function
    y = myf;
end

function [c,ceq] = constr(x,Nom_values,new_corr,K,theta1,theta2,epsi1,epsi2,w,pci_func,pci_den,para_selection,sampling,bio_nom,sf_nom_mod)
    if ~isequal(x,xLast) % Check if computation is necessary
        [myf,myc,myceq] = compute_robust_loop_input_uncertainty(x,Nom_values,new_corr,K,theta1,theta2,epsi1,epsi2,w,pci_func,pci_den,para_selection,sampling,bio_nom,sf_nom_mod);
        xLast = x;
    end
    % Now compute constraint functions
    c = myc; % In this case, the computation is trivial
    ceq = myceq;
end

disp([-myf myc+120])
end