function [U_doe,grad_dir_out,grad_step_out,lb,ub] = experimental_design_pen(K_model,initial_cond,U_vec,sampling,extra_var,...
    U_out_vec,i_b2b,num_grad_selec,para_selection,num_grad_points,grad_point_vec,...
    grad_dir_vec,grad_step_size_vec,pen_std)
%function to compute the next batch run used for gradient measurements. 
%Author: Rubin Hille


%define model
model = @penicillin_process_model;

Y0 = initial_cond;
Y0_grad1 = Y0;
Y0_grad2 = Y0;
U_curr = U_vec(end,:);

U = U_vec(i_b2b,:);

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

num_para = length(K_model);

K_ini = K_model;

%calculate gradient sensitivities from output measurements from past
%operating points

estimation_crit = 'simulation';
time_vec = sampling;
vcd_vec = [];

n_para_doe = 4;


        
        plant_grad_vec = [];
        grad_pred_vec = [];
        
%         K_new = K_ini;
        %only change the parameters which have been selected for gradient
        %correction
        
%         for i_selection = 1:length(K_ini)
%           if para_selection(i_selection) == 1
%               K_new(i_selection) = x(i_selection);
%           end
%         end
        
        %calculate for nominal conditions
%         U_curr = U_out_vec(end,:);
        U(1) = U_curr(1);
        Y0_grad1(3) = U(1);
        %get experimental data of corresponding batch
%         proc_out_mean = reshape(proc_out_mean_history(:,i_b2b),length(sampling),num_outputs);
        
        %get vcd data of corresponding batch (for interpolation)
        time_vec = [];
        vcd_vec = [];
        [~, Y_nom] = ode45(model,sampling, Y0_grad1, opt,U,K_ini,extra_var);
    
        %nominal cost
        Obj_curr = Y_nom(end,2)*Y_nom(end,4);
        
        %using the last "num_corrected_grad" gradient measurements 
        i_1 = 1;
        for i_corr_grad = 1:5
%             i_corr_grad = 1:length(grad_point_vec)
            
            if grad_point_vec(end-i_corr_grad+1) == 1
                %get deviation between correspodning operating point and 
                U_new = U_out_vec(end-i_corr_grad,:);
                U_dev_new = U_new-U_curr;
                
                %calculate predicted gradient
                U(1) = U_new(1);
                Y0_grad2(3) = U(1);
                [~, Y_grad] = ode45(model,sampling, Y0_grad2, opt,U,K_ini,extra_var);

                %compute nominal gradient
                grad_pred_nominal = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/sqrt(sum(U_dev_new.^2)); 
                
                %calculate gradients with a deviation in parameter values
                i_2 = 1;
                for i_para = 1:num_para

                 %only calculate sensitivity for parameters which are used for
                 %gradient correction
                  if para_selection(i_para) == 1

                    %define new parameter values
                    K_new = K_ini;

                    %deviation in parameter for sensitivity step size
                    dev_para = sqrt(eps)*(K_new(i_para)+0.1);

                    %change para value by small deviation
                    K_new(i_para) = K_new(i_para)+dev_para;

                    %calculate for nominal conditions
           
%                     Y0_grad1(4) = U_curr(1);
%                     [~, Y_nom_new] = ode45(model,sampling, Y0_grad1, opt,K_new,estimation_crit,time_vec,vcd_vec);
                    U(1) = U_curr(1);
                    Y0_grad1(3) = U(1);
                    [~, Y_nom_new] = ode45(model,sampling, Y0_grad1, opt,U,K_new,extra_var);
                    
                    U(1) = U_new(1);
                    Y0_grad2(3) = U(1);
                    [~, Y_grad_new] = ode45(model,sampling, Y0_grad2, opt,U,K_new,extra_var);
                    
%                     Y0_grad2(4) = U_new(1);
%                     [~, Y_grad_new] = ode45(model,sampling, Y0_grad2, opt,K_new,estimation_crit,time_vec,vcd_vec);

                    %compute gradient
                    grad_pred_new = (Y_grad_new(end,2)*Y_grad_new(end,4)-Y_nom_new(end,2)*Y_nom_new(end,4))/sqrt(sum(U_dev_new.^2));
                    
                    %calculate gradient sensitivities
                    grad_sens = (grad_pred_new-grad_pred_nominal)./abs(dev_para);

                    %scale with nominal parameter and gradient value
%                     grad_sens =  grad_sens./abs(grad_pred_nominal);
                    grad_sens =  grad_sens.*K_new(i_para);

                    %put gradient sensitivities into matrix
                    grad_sens_matrix(i_1,i_2) = grad_sens;

                    i_2 = i_2 + 1;
                    
                    %only use small number of parameters
                    if i_2 >= n_para_doe
                        break;
                    end
                    
                  end
    
                end
            
            covar_vec(i_1) = norm(U_dev_new)^2;    
            i_1 = i_1 + 1;
            end
            
            %test with only num_grad_points-1 operating points
            if i_1 == num_grad_points-2
%               i_1 == num_grad_points-1
                break;
            end
         
        end
        
        
%get the initial gradient sensitivity matrix
gradient_sens_matrix_initial = grad_sens_matrix;

%gradient covariance matrix initial
initial_covar_mat = diag(covar_vec);

%calculate gradient sensitivities for a range of operating points

%range of inputs
U_new_grad_points = linspace(1,75,30);

        U(1) = U_curr(1);
        Y0_grad1(3) = U(1);
        %get experimental data of corresponding batch
%         proc_out_mean = reshape(proc_out_mean_history(:,i_b2b),length(sampling),num_outputs);
        
        %get vcd data of corresponding batch (for interpolation)
        time_vec = [];
        vcd_vec = [];
        [~, Y_nom] = ode45(model,sampling, Y0_grad1, opt,U,K_ini,extra_var);
        
        %using the last "num_corrected_grad" gradient measurements 
        i_1 = 1;
        for i_corr_grad = 1:length(U_new_grad_points)
            
                %get deviation between correspodning operating point and 
                U_new = U_new_grad_points(i_corr_grad);
                U_dev_new = U_new-U_curr;
                
                %calculate predicted gradient
                U(1) = U_new(1);
                Y0_grad2(3) = U(1);
                [~, Y_grad] = ode45(model,sampling, Y0_grad2, opt,U,K_ini,extra_var);
%                 Y0_grad2(4) = U_new(1);
%                 [~, Y_grad] = ode45(model,sampling, Y0_grad2, opt,K_ini,estimation_crit,time_vec,vcd_vec);

                %get objective function value at U_new
                Obj_function_value(i_corr_grad) = Y_grad(end,2)*Y_grad(end,4);

                %compute nominal gradient
%                 grad_pred_nominal = (Y_grad(end,11)-Y_nom(end,11))/sqrt(sum(U_dev_new.^2)); 
                grad_pred_nominal = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/sqrt(sum(U_dev_new.^2)); 
                
                %calculate gradients with a deviation in parameter values
                i_2 = 1;
                for i_para = 1:num_para

                 %only calculate sensitivity for parameters which are used for
                 %gradient correction
                  if para_selection(i_para) == 1

                    %define new parameter values
                    K_new = K_ini;

                    %deviation in parameter for sensitivity step size
                    dev_para = sqrt(eps)*(K_new(i_para)+0.1);

                    %change para value by small deviation
                    K_new(i_para) = K_new(i_para)+dev_para;

                    %calculate for nominal conditions
           
%                     Y0_grad1(4) = U_curr(1);
%                     [~, Y_nom_new] = ode45(model,sampling, Y0_grad1, opt,K_new,estimation_crit,time_vec,vcd_vec);
%                     
%                     Y0_grad2(4) = U_new(1);
%                     [~, Y_grad_new] = ode45(model,sampling, Y0_grad2, opt,K_new,estimation_crit,time_vec,vcd_vec);
                    
                    U(1) = U_curr(1);
                    Y0_grad1(3) = U(1);
                    [~, Y_nom_new] = ode45(model,sampling, Y0_grad1, opt,U,K_new,extra_var);
                    
                    U(1) = U_new(1);
                    Y0_grad2(3) = U(1);
                    [~, Y_grad_new] = ode45(model,sampling, Y0_grad2, opt,U,K_new,extra_var);

                    %compute gradient
                    grad_pred_new = (Y_grad_new(end,2)*Y_grad_new(end,4)-Y_nom_new(end,2)*Y_nom_new(end,4))/sqrt(sum(U_dev_new.^2));

%                     grad_pred_new = (Y_grad_new(end,11)-Y_nom_new(end,11))/sqrt(sum(U_dev_new.^2)); 
                    
                    %calculate gradient sensitivities
                    grad_sens = (grad_pred_new-grad_pred_nominal)./abs(dev_para);

                    %scale with nominal parameter and gradient value
%                     grad_sens =  grad_sens./abs(grad_pred_nominal);
                    grad_sens =  grad_sens.*K_new(i_para);

                    %put gradient sensitivities into matrix
                    grad_sens_matrix2(:,i_2) = grad_sens;

                    i_2 = i_2 + 1;
                    
                    %only use small number of parameters
                    if i_2 >= n_para_doe
                        break;
                    end
                    
                  end
                
                
                end
        
        %initial covariance vector        
        covar_vec_new = covar_vec;
        
        %new element
        covar_new = norm(U_dev_new)^2;  
        
        %update covariance vector
        covar_vec_new = [covar_vec_new covar_new];
        
        %covar matrix
        covar_mat_new = diag(covar_vec_new);
        
        %gradient sensitivity from new point   
        new_grad_sens = grad_sens_matrix2;
        
        %compute new gradient sensitivity matrix
        new_gradient_matrix = [gradient_sens_matrix_initial; new_grad_sens];    
%         new_gradient_matrix = [gradient_sens_matrix_initial(end,:); new_grad_sens];  
%         new_gradient_matrix = [new_grad_sens]; 
        
        %calculate FIM
        FIM_new = new_gradient_matrix'*inv(covar_mat_new)*new_gradient_matrix;
        
        %D optimality criterion
        DoE_crit(i_corr_grad) = det(FIM_new+0*eye(size(FIM_new)));
%         /(1+0.01*abs(U_new(1)-U_curr(1)))
        
        %adjusted for distance
%         alpha = -(mean(abs(DoE_crit(1))))/(norm(U_new_grad_points(1) - U_curr))/2;
%         DoE_crit_2(i_corr_grad) = DoE_crit(i_corr_grad)+alpha*norm(U_new - U_curr);

        %modified E criterion
%         E = eig(inv(FIM_new+0*eye(size(FIM_new))));
%         DoE_crit(i_1) = 1/(max(E)/min(E));
         
%          i_1 = i_1 + 1;
        end

%adjust DOE criterion for values close to the current operating point
% DoE_crit = DoE_crit.*(abs(U_new_grad_points - U_curr) > 3);

%get DoE_crit at U_curr
% DoE_U_Curr = interp1(U_new_grad_points,DoE_crit,U_curr);
% alpha = -abs(DoE_U_Curr)/norm(U_new_grad_points(1)-U_curr)/4;
% for i_doe_2 = 1:length(DoE_crit)
%     DoE_crit_2(i_doe_2) = DoE_crit(i_doe_2)+alpha*norm(U_new_grad_points(i_doe_2)-U_curr);
% end

figure(15)
clf
scatter(U_new_grad_points,DoE_crit,'b','fill');
% hold on
% scatter(U_new_grad_points,DoE_crit_2,'r','fill');
grid on;
drawnow;

%% Determine range where loss is less than 0.1 for example
Obj_function_value_rel = Obj_function_value/Obj_curr;

figure(17)
clf
scatter(U_new_grad_points,Obj_function_value_rel,'k','fill')
grid on;
drawnow;

if i_b2b > 2
    %cost eps from standard deviation
    cost_eps = 1 - pen_std/Obj_curr
else
    cost_eps = 0.90;
end


%go through vectors to determine lower and upper bound
new_grad_points2 = linspace(1,105,65);
new_obj_value = interp1(U_new_grad_points,Obj_function_value_rel,new_grad_points2);


cost_bound_counter = 1;
for i_cost_bound = 1:length(new_obj_value)
    if cost_bound_counter == 1 && new_obj_value(i_cost_bound) >= cost_eps
        lb = new_grad_points2(i_cost_bound);
        cost_bound_counter = 2;
        %in case no upper bound found
        ub = lb + 15;
    elseif cost_bound_counter == 2 && new_obj_value(i_cost_bound) <= cost_eps
        ub = new_grad_points2(i_cost_bound);
        cost_bound_counter = 3;
    end
end

% lb
ub = ub + 2.7;


%% find most informative point within range

% % if i_b2b > 1
%limit on permissible distance from current operating point
delta_u_limit = 15;

U0 = U_curr+grad_step_size_vec(1);

% % lb = max(0.1,U_curr(1) - delta_u_limit);
% % ub = U_curr(1) + delta_u_limit;

% options=optimset('Display','iter',...
%     'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'Algorithm','interior-point');
% [x,f,~,~] = fmincon(@DoE_obj,U0,[],[],[],[],lb,ub,[],options);


%evaluate doe objective at lb and ub
F_lb = DoE_obj(lb);
F_ub = DoE_obj(ub);

if F_lb < F_ub
    x = lb;
else
    x = ub;
end

% %test with fminbnd
% options=optimset('Display','iter');
% [x,f,~,~] = fminbnd(@DoE_obj,lb,ub,options);
% 
% 
% if abs(x-U_curr(1)) < 0.5
%     x = U_curr(1) + sign(x)*0.5;
% end
% 
% 
% if sum(abs(x-U_out_vec) < 0.5) >= 1
%     x = x + sign(U_curr(1)-x)*1;
% end

disp('CURRENT OPERATING POINT')
disp(U_curr(1))
disp('NEW EXPERIMENTAL GRADIENT POINT')
U_doe = x;
disp(U_doe)

grad_dir_out = 1;
grad_step_out = U_doe-U_curr(1);


% % else
% %     
% % U_doe = 0;
% % grad_dir_out=0;
% % grad_step_out=0;
% % end

    function F = DoE_obj(x)
        
        U(1) = U_curr(1);
        Y0_grad1(3) = U(1);
        %get experimental data of corresponding batch
%         proc_out_mean = reshape(proc_out_mean_history(:,i_b2b),length(sampling),num_outputs);
        
        %get vcd data of corresponding batch (for interpolation)
        time_vec = [];
        vcd_vec = [];
        [~, Y_nom] = ode45(model,sampling, Y0_grad1, opt,U,K_ini,extra_var);
        
        %using the last "num_corrected_grad" gradient measurements 
        i_1 = 1;
        for i_corr_grad = 1:1
            
                %get deviation between correspodning operating point and 
                U_new = x;
                U_dev_new = U_new-U_curr;
                
                %calculate predicted gradient
                U(1) = U_new(1);
                Y0_grad2(3) = U(1);
                [~, Y_grad] = ode45(model,sampling, Y0_grad2, opt,U,K_ini,extra_var);

                %compute nominal gradient
                grad_pred_nominal = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/sqrt(sum(U_dev_new.^2)); 

                
                %calculate gradients with a deviation in parameter values
                i_2 = 1;
                for i_para = 1:num_para

                 %only calculate sensitivity for parameters which are used for
                 %gradient correction
                  if para_selection(i_para) == 1

                    %define new parameter values
                    K_new = K_ini;

                    %deviation in parameter for sensitivity step size
                    dev_para = sqrt(eps)*(K_new(i_para)+0.1);

                    %change para value by small deviation
                    K_new(i_para) = K_new(i_para)+dev_para;

                    %calculate for nominal conditions
                    U(1) = U_curr(1);
                    Y0_grad1(3) = U(1);
                    [~, Y_nom_new] = ode45(model,sampling, Y0_grad1, opt,U,K_new,extra_var);
                    
                    U(1) = U_new(1);
                    Y0_grad2(3) = U(1);
                    [~, Y_grad_new] = ode45(model,sampling, Y0_grad2, opt,U,K_new,extra_var);

                    %compute gradient
                    grad_pred_new = (Y_grad_new(end,2)*Y_grad_new(end,4)-Y_nom_new(end,2)*Y_nom_new(end,4))/sqrt(sum(U_dev_new.^2));
                    
                    %calculate gradient sensitivities
                    grad_sens = (grad_pred_new-grad_pred_nominal)./abs(dev_para);

                    %scale with nominal parameter and gradient value
%                     grad_sens =  grad_sens./abs(grad_pred_nominal);
                    grad_sens =  grad_sens.*K_new(i_para);

                    %put gradient sensitivities into matrix
                    grad_sens_matrix3(:,i_2) = grad_sens;

                    i_2 = i_2 + 1;
                    
                    %only use small number of parameters
                    if i_2 >= n_para_doe
                        break;
                    end
                    
                  end
                
                
                end
        end
        
        %initial covariance vector        
        covar_vec_new = covar_vec;
        
        %new element
        covar_new = norm(U_dev_new)^2;  
        
        %update covariance vector
        covar_vec_new = [covar_vec_new covar_new];
        
        %covar matrix
        covar_mat_new = diag(covar_vec_new);
        
        %gradient sensitivity from new point        
        new_grad_sens = grad_sens_matrix3;    
        
        %compute new gradient sensitivity matrix
        new_gradient_matrix = [gradient_sens_matrix_initial; new_grad_sens];    
%         new_gradient_matrix = [gradient_sens_matrix_initial(end,:); new_grad_sens];  
%         new_gradient_matrix = [new_grad_sens]; 
        
        %calculate FIM
        FIM_new = new_gradient_matrix'*inv(covar_mat_new)*new_gradient_matrix;
        
        %D optimality criterion
        DoE_crit = det(FIM_new+0*eye(size(FIM_new)));
%         /(1+0.01*abs(U_new(1)-U_curr(1)))
        
        F = -DoE_crit;
        
    end

end