function [process_optimum] = optimize_penicillin_process(initial_cond,...
    initial_decision_val,sampling_time,extra_var,vol_constr,K_process)

%function to optimize the "true" penicillin process


%define upper and lower bounds
lb=[0.1 0.01];
ub=[100 100];

%intial values for initial substrate concentration and flowrate
U0 = initial_decision_val;

%define process simulator
pen_simulator = @penicillin_process_simulator;


%sampling time
sampling= sampling_time;

%define feed concentration
sf_nom = extra_var(1);

%constraint on the reactor volum
reactor_volume_constr = vol_constr;

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% process optimization

disp('PENICILLIN "TRUE" PROCESS OPTIMIZATION');
options=optimset('Display','iter',...
    'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'Algorithm','active-set');
[x,f,~,~] = fmincon(@objfun,U0,[],[],[],[],lb,ub,@constr,options);

process_optimum = x;

disp('Optimal Initial Substrate Concentration');
disp(x(1));
disp('Optimal Feed Rate');
disp(x(2));

%objective function for maximizing the amount of penicillin at the end
%of the batch
function y = objfun(x)
    
    [myf,myc,myceq] = sim_process(x);
    
    y = myf;
end

function [c,ceq] = constr(x)
   
    
    [myf,myc,myceq] = sim_process(x);
    % Now compute constraint functions
    c = myc; % In this case, the computation is trivial
    ceq = myceq;
end

    function [myf,myc,myceq] = sim_process(x1)
        
    U = x1;
    
    %initial process conditions
    Y0 = initial_cond;
    %update initial substrate concentration
    Y0(3) = U(1);
    

    [~, Y] = ode15s(pen_simulator,[sampling], Y0, opt,U,K_process,sf_nom);
%     
%     S_inter = interp1(sampling,Y(:,2),U(3))*interp1(sampling,Y(:,4),U(3));
%     constr_inter = interp1(sampling,Y(:,4),U(3));
    
    S=Y(end,2)*Y(end,4);
    constr=Y(end,4);

% 
%     S = S_inter;
%     constr = constr_inter;

    f1=-S;
    c1=constr-reactor_volume_constr;
    ceq1=[];
    
    myf = f1;
    myc = c1;
    myceq = ceq1;
    end
end