function [process_opt] = optimize_simple_process(initial_decision_val,sampling,...
    Y_initial,Para_sim,cost_para)


%define upper and lower bounds
lb = 0;
ub = 100;

%intial values for initial concentration
U = initial_decision_val;

%define model used for process simulation
simple_simulator = @simple_process_simulator;


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
process_opt = x;


    function F = objfun(x)
        
        %run numerical ode solver 
        [T, Y] = ode45(simple_simulator,sampling,Y0,opt,Para_sim,x);
        
        %cost function
        F = -(Y(end) - c*x^2);
        
    end

end