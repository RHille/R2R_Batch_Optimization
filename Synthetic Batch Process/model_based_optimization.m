function [U_new] = model_based_optimization(U,sampling,...
    Y_initial,para_model,cost_para)

%define upper and lower bounds
lb = 0;
ub = 100;


%define model used for process simulation
simple_model = @simple_process_model;


opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% process optimization

Y0 = Y_initial;

%parameter of cost function
c = cost_para;

disp('SIMPLE PROCESS OPTIMIZATION');
options=optimset('Display','iter',...
    'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'Algorithm','active-set');
[x,f,~,~] = fmincon(@objfun,U,[],[],[],[],lb,ub,[],options);

%process optimum
U_new = x;


    function F = objfun(x)
        
        %run numerical ode solver 
        [T, Y] = ode45(simple_model,sampling,Y0,opt,para_model(1:2),x);
        
        %cost function
        F = -(Y(end) - para_model(3)*x^2);
        
    end


end