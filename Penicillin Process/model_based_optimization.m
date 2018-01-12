function [U_new] = model_based_optimization(U_vec,K_model_adapt,case_study,...
    lb,ub,process_constr,extra_var,initial_cond,new_corr,sampling,i_b2b)
%function to perform model based optimization

disp('MODEL-BASED OPTIMIZATION');

%define model
if strcmp(case_study,'penicillin')
    %Define penicillin process simulator
    model = @penicillin_process_model;
end



%initial process conditions
Y0 = initial_cond;

%ode options
opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% 'NonNegative',[1,2,3,4],

U_ini = U_vec(i_b2b,:);


%fmincon options
options=optimset('Algorithm','interior-point','Display','iter','UseParallel','always',...
    'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6);

%'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-8
[x,fval] = fmincon(@objfun,U_ini,[],[],[],[],lb,ub,@constrfun,options);

U_new = x;

function y = objfun(x)
    
    %call up model
    [myf,myc,myceq] = sim_model(x);
    
    y = myf;
end

function [c,ceq] = constrfun(x)
   
    %call up model
    [myf,myc,myceq] = sim_model(x);
    
    % Now compute constraint functions
    c = myc; 
    ceq = myceq;
end

    function [myf,myc,myceq] = sim_model(x1)
        
    U = x1;
    
    %update initial substrate concentration
    Y0(3) = U(1);
    

    [~, Y] = ode45(model,[sampling], Y0, opt,U,K_model_adapt,extra_var);
%     
%     S_inter = interp1(sampling,Y(:,2),U(3))*interp1(sampling,Y(:,4),U(3));
%     constr_inter = interp1(sampling,Y(:,4),U(3));

    Y = Y - new_corr;
    
    %Final amount of penicillin at the end of batch
    S=Y(end,2)*Y(end,4);
    
    %Volume at the end of the batch
    constr=Y(end,4);

% 
%     S = S_inter;
%     constr = constr_inter;

    f1=-S;
    c1=constr-process_constr;
    ceq1=[];
    
    myf = f1;
    myc = c1;
    myceq = ceq1;
    end

end