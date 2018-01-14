function [para_corr,new_corr,gradient_model,gradient_new] = gradient_correction(U,sampling,...
    Y_initial,para_model,gradient,prev_corr,perturbation,c_process,c_model,t_error,...
    U_out_vec,Obj_fun_vec,num_grad_corr,i_b2b,grad_point_vec)


% %define upper and lower bounds
% lb = 0.1*ones(size(para_model));
% ub = 5*ones(size(para_model));

delta = 0.5;
lb = para_model-delta*abs(para_model);
ub = para_model+delta*abs(para_model);

num_ref_points = 0;

%define model used for process simulation
simple_model = @simple_process_model;


opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% process optimization

Y0 = Y_initial;

%initial parameter values
K_ini = para_model;


%simulate process
[T, y_model_old] = ode45(simple_model,sampling,Y0,opt,para_model(1:2),U);

%calculate sensitivities
for i_sens = 1:length(K_ini)
            
        %update parameter value
        para_new = K_ini;

        para_dev = sqrt(eps)*(para_new(i_sens));
        para_new(i_sens) = para_new(i_sens) + para_dev;

        [~, y_dev] = ode45(simple_model,sampling,Y0,opt,para_new,U);

        sens(i_sens,:) = (y_dev-y_model_old)/para_dev;
end


disp('GRADIENT CORRECTION');
%fmincon options
options=optimset('Algorithm','active-set','Display','iter','UseParallel','always',...
    'TolX',1e-4,'TolFun',1e-4,'TolCon',1e-4,'MaxFunEvals',1000);

%'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-8
[x,SSE] = fmincon(@gradobj_ext,K_ini,[],[],[],[],lb,ub,@gradconstr,options);

% @gradconstr

% x(3) = K_ini(3);

%new parameter values
para_corr = x;

%calculate correction term
[T, y_model_new] = ode45(simple_model,sampling,Y0,opt,para_corr,U);
c_new = -y_model_old+y_model_new;
new_corr=prev_corr+c_new;
    

    function F = gradobj_ext(x)
        
        
%         x(3) = K_ini(3);
       
        for i2 = 0:min(length(U_out_vec)-1,num_ref_points)

        U_curr = U_out_vec(end);
        Obj_curr = Obj_fun_vec(end);

        [T, y_model_nominal] = ode45(simple_model,sampling,Y0,opt,x(1:2),U_curr);
        
        
        for i=1:length(U_out_vec)-1
            
%             min(num_grad_corr,length(U_out_vec)-1)
            
            if grad_point_vec(end-i) == 1
            
            %get inputs from vector
            U_new = U_out_vec(end-i);
            Obj_new = Obj_fun_vec(end-i);
            dev = abs(U_curr-U_new);
            gradient_new = (Obj_new-Obj_curr)/(max(1e-3,dev));
            
                grad_error(i) = gradobj(x,U_new,y_model_nominal,dev,gradient_new,U_curr);            
            else 
                grad_error(i) = 0;
            end
        end

            F1(i2+1) = sum(grad_error.^2);
        end
        
        F = sum(F1);
        
    end
    
    function F = gradobj(x,U_new,y_model_nominal,dev,gradient_new,U_curr)
        
        %calculate predicted gradient
        
        
%         %nominal
%         [T, y_model_nominal] = ode45(simple_model,sampling,Y0,opt,x,U);
        
%         U_perturbed = U + perturbation;
        [T, y_model_perturbed] = ode45(simple_model,sampling,Y0,opt,x(1:2),U_new);
        

        %calculate nominal cost
        cost_nominal = y_model_nominal(end) - x(3)*U_curr^2;

        %perturbed cost
        cost_perturbed = y_model_perturbed(end) - x(3)*U_new^2;

        %estimate gradient
        gradient_model = (cost_perturbed - cost_nominal)/max(1e-3,dev);
        
        %minimize the difference between measured and predicted gradient

        F = norm(gradient_new - gradient_model);
%         
        
    end

    function [myc,myceq] = gradconstr(x)
        
%         x(3) = K_ini(3);
        
        para_nom = x;
        [~, y_constr] = ode45(simple_model,sampling,Y0,opt,para_nom(1:2),U);

        c = 0;
        for i_match = 1:length(x)
                c = c + sens(i_sens,:)*(para_nom(i_match)-K_ini(i_match));
        end
        
        Y_constr_prime=y_constr-c'-prev_corr;  %new model
        Y_prime=y_model_old-prev_corr;  %previous model
        trunc_e=abs(Y_constr_prime./Y_prime-1);
        
        myc(1)= max(max(trunc_e(:,:)))-t_error;
        
%         myc(1)= max(trunc_e(end))-t_error;
        
        myceq = [];
    end


end