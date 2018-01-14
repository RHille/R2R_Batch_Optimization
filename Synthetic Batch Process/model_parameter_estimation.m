function [para_new] = model_parameter_estimation(U,sampling,...
    Y_initial,para_model,y_out,prev_corr,y_max,y_min)



%define upper and lower bounds
lb = 0.1*ones(size(para_model));
ub = 5*ones(size(para_model));

% delta = 0.25;
% lb = para_model-delta*abs(para_model);
% ub = para_model+delta*abs(para_model);


%define model used for process simulation
simple_model = @simple_process_model;


opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% process optimization

Y0 = Y_initial;

%initial parameter values
K_ini = para_model;


%fmincon options
options=optimset('Algorithm','active-set','Display','iter','UseParallel','always',...
    'TolX',1e-4,'TolFun',1e-4,'TolCon',1e-4);

%'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-8
[x,SSE] = fmincon(@sseobj,K_ini,[],[],[],[],lb,ub,@sb_constr,options);

%@sb_constr
%new parameter values
x(3) = K_ini(3);
para_new = x;

    function F = sseobj(x)
        
        x(3) = K_ini(3);
        %get model prediction
        [T, y_model] = ode45(simple_model,sampling,Y0,opt,x(1:2),U);
        y_model = y_model - prev_corr;
        
        %calculate sse objective
        sse = sum((y_out-y_model).^2);
%         +0*1e1*norm(x(1:2)-K_ini)^2
        F = sse ;
%         + 1e3*norm(x(1:2)-K_ini(1:2))^2
    end

    function [myc,myceq] = sb_constr(x)
        
         x(3) = K_ini(3);
        %run the "model" with initial values     
        [T, y_model] = ode45(simple_model,sampling,Y0,opt,x(1:2),U);
        %add correction for accurate prediction
        y_model = y_model - prev_corr;
           
        myc = [y_model - y_max';y_min' - y_model];
        
        myceq = 0;
    end


end