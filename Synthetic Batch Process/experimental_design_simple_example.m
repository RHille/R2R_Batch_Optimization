function [doe_out,lb_doe,ub_doe] = experimental_design_simple_example(U,sampling,...
    Y_initial,para_model,grad_point_vec,U_out_vec,num_grad_points,cost_std)


disp('EXPERIMENTAL DESIGN');

%define model
simple_model = @simple_process_model;

Y0 = Y_initial;
Y0_grad1 = Y0;
Y0_grad2 = Y0;
U_curr = U;

K_model = para_model;

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

num_para = length(K_model);

K_ini = K_model;

num_corrected_grad_points = min(5,num_grad_points);


para_selection = [1 1 1];

%calculate gradient sensitivities from output measurements from past
%operating points
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
        U = U_curr;
        [~, Y_nom] = ode45(simple_model,sampling,Y0_grad1,opt,K_ini(1:2),U);
    
        %nominal cost
        Obj_curr = Y_nom(end) - K_ini(3)*U^2;
        
        %using the last "num_corrected_grad" gradient measurements 
        i_1 = 1;
        for i_corr_grad = 1:num_corrected_grad_points
%             i_corr_grad = 1:length(grad_point_vec)

            if grad_point_vec(end-(i_corr_grad-1)) == 1
                %get deviation between correspodning operating point and 
                U_new = U_out_vec(end-(i_corr_grad-1),:);
                U_dev_new = U_new-U_curr;
                
                %calculate predicted gradient
                U(1) = U_new(1);
                [~, Y_grad] = ode45(simple_model,sampling, Y0_grad2,opt,K_ini(1:2),U);
  
                Obj_new = Y_grad(end) - K_ini(3)*U^2;
                %compute nominal gradient
                grad_pred_nominal = (Obj_new-Obj_curr)/sqrt(sum(U_dev_new.^2)); 
                
                %calculate gradients with a deviation in parameter values
                i_2 = 1;
                for i_para = 1:num_para

                 %only calculate sensitivity for parameters which are used for
                 %gradient correction
                  if para_selection(i_para) == 1

                    %define new parameter values
                    K_new = K_ini;

                    %deviation in parameter for sensitivity step size
                    dev_para = sqrt(eps)*(K_new(i_para)+0.01);

                    %change para value by small deviation
                    K_new(i_para) = K_new(i_para)+dev_para;

                    %calculate for nominal conditions
           
%                     Y0_grad1(4) = U_curr(1);
%                     [~, Y_nom_new] = ode45(model,sampling, Y0_grad1, opt,K_new,estimation_crit,time_vec,vcd_vec);
                    U(1) = U_curr(1);
                    [~, Y_nom_new] = ode45(simple_model,sampling,Y0_grad1,opt,K_new(1:2),U);
                    
                    Obj_curr_new = Y_nom_new(end) - K_new(3)*U^2;
                    
                    U(1) = U_new(1);
                    [~, Y_grad_new] = ode45(simple_model,sampling, Y0_grad2,opt,K_new(1:2),U);
                    
                    Obj_new_new = Y_grad_new(end) - K_new(3)*U^2;
                    
%                     Y0_grad2(4) = U_new(1);
%                     [~, Y_grad_new] = ode45(model,sampling, Y0_grad2, opt,K_new,estimation_crit,time_vec,vcd_vec);

                    %compute gradient
                    grad_pred_new = (Obj_new_new-Obj_curr_new)/sqrt(sum(U_dev_new.^2)); 
                    
                    %calculate gradient sensitivities
                    grad_sens = (grad_pred_new-grad_pred_nominal)./abs(dev_para);

                    %scale with nominal parameter and gradient value
%                     grad_sens =  grad_sens./abs(grad_pred_nominal);
                    grad_sens =  grad_sens.*K_new(i_para);

                    %put gradient sensitivities into matrix
                    grad_sens_matrix(i_1,i_2) = grad_sens;

                    i_2 = i_2 + 1;
                    
% %                     %only use small number of parameters
% %                     if i_2 >= n_para_doe
% %                         break;
% %                     end
                    
                  end
    
                end
            
%             covar_vec(i_1) = 1/norm(U_dev_new)^2;   
            covar_vec(i_1) = 1; 
            i_1 = i_1 + 1;
            end
            
%             %test with only num_grad_points-1 operating points
%             if i_1 == num_grad_points-2
% %               i_1 == num_grad_points-1
%                 break;
%             end
         
        end
        
        
%get the initial gradient sensitivity matrix
gradient_sens_matrix_initial = grad_sens_matrix;

%gradient covariance matrix initial
initial_covar_mat = diag(covar_vec);

%calculate gradient sensitivities for a range of operating points

%range of inputs
U_new_grad_points = linspace(1,20,40);

        U(1) = U_curr(1);
        %get experimental data of corresponding batch
%         proc_out_mean = reshape(proc_out_mean_history(:,i_b2b),length(sampling),num_outputs);
        
        %get vcd data of corresponding batch (for interpolation)
        [~, Y_nom] = ode45(simple_model,sampling,Y0_grad1,opt,K_ini(1:2),U);
        
        Obj_function_value =  Y_nom(end) - K_ini(3)*U^2;
        obj_fun_value_curr = Obj_function_value;
        
        %using the last "num_corrected_grad" gradient measurements 
        i_1 = 1;
        for i_corr_grad = 1:length(U_new_grad_points)
            
                %get deviation between correspodning operating point and 
                U_new = U_new_grad_points(i_corr_grad);
                U_dev_new = U_new-U_curr;
                
                %calculate predicted gradient
                U(1) = U_new(1);
                [~, Y_grad] = ode45(simple_model,sampling, Y0_grad2, opt,K_ini(1:2),U);
%                 Y0_grad2(4) = U_new(1);
%                 [~, Y_grad] = ode45(model,sampling, Y0_grad2, opt,K_ini,estimation_crit,time_vec,vcd_vec);

                %get objective function value at U_new
                Obj_function_value_grad =  Y_grad(end) - K_ini(3)*U^2;
                obj_fun_value_vec(i_corr_grad) = Obj_function_value_grad;

                %compute nominal gradient
%                 grad_pred_nominal = (Y_grad(end,11)-Y_nom(end,11))/sqrt(sum(U_dev_new.^2)); 
                grad_pred_nominal = (Obj_function_value_grad-Obj_function_value)/sqrt(sum(U_dev_new.^2));
                
                %calculate gradients with a deviation in parameter values
                i_2 = 1;
                grad_sens_matrix2 = [];
                for i_para = 1:num_para

                 %only calculate sensitivity for parameters which are used for
                 %gradient correction
                  if para_selection(i_para) == 1

                    %define new parameter values
                    K_new = K_ini;

                    %deviation in parameter for sensitivity step size
                    dev_para = sqrt(eps)*(K_new(i_para)+0.01);

                    %change para value by small deviation
                    K_new(i_para) = K_new(i_para)+dev_para;

                    %calculate for nominal conditions
           
%                     Y0_grad1(4) = U_curr(1);
%                     [~, Y_nom_new] = ode45(model,sampling, Y0_grad1, opt,K_new,estimation_crit,time_vec,vcd_vec);
%                     
%                     Y0_grad2(4) = U_new(1);
%                     [~, Y_grad_new] = ode45(model,sampling, Y0_grad2, opt,K_new,estimation_crit,time_vec,vcd_vec);
                    
                    U(1) = U_curr(1);
%                     Y0_grad1(3) = U(1);
                    [~, Y_nom_new] = ode45(simple_model,sampling, Y0_grad1, opt,K_new(1:2),U);
                   
                    Obj_function_value_new =  Y_nom_new(end) - K_new(3)*U^2;
                    U(1) = U_new(1);
%                     Y0_grad2(3) = U(1);
                    [~, Y_grad_new] = ode45(simple_model,sampling, Y0_grad2, opt,K_new(1:2),U);
                    
                    Obj_function_value_grad_new =  Y_grad_new(end) - K_new(3)*U^2;

                    %compute gradient
                    grad_pred_new = (Obj_function_value_grad_new-Obj_function_value_new)/sqrt(sum(U_dev_new.^2));

%                     grad_pred_new = (Y_grad_new(end,11)-Y_nom_new(end,11))/sqrt(sum(U_dev_new.^2)); 
                    
                    %calculate gradient sensitivities
                    grad_sens = (grad_pred_new-grad_pred_nominal)./abs(dev_para);

                    %scale with nominal parameter and gradient value
%                     grad_sens =  grad_sens./abs(grad_pred_nominal);
                    grad_sens =  grad_sens.*K_new(i_para);

                    %put gradient sensitivities into matrix
%                     grad_sens_matrix2(:,i_2) = grad_sens
                  
                    grad_sens_matrix2(i_2) = grad_sens;

                    i_2 = i_2 + 1;
                    
                    %only use small number of parameters
% %                     if i_2 >= n_para_doe
% %                         break;
% %                     end
                    
                  end
                
                
                end
        
        %initial covariance vector        
        covar_vec_new = covar_vec;
        
        %new element
        covar_new = 1/norm(U_dev_new)^2; 
        covar_new = 1; 
        
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
        DoE_crit(i_corr_grad) = det(FIM_new+0.1*eye(size(FIM_new)));
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

% figure(15)
% clf
% scatter(U_new_grad_points,DoE_crit,'b','fill');
% % hold on
% % scatter(U_new_grad_points,DoE_crit_2,'r','fill');
% grid on;
% drawnow;

%calculate minimum distance to reduce gradient uncertainty
delta_obj = obj_fun_value_curr - obj_fun_value_vec;
% cost_std
% [U_new_grad_points' delta_obj']
delta_count = 1;
for i2 = 1:length(delta_obj)
    if abs(delta_obj(i2)) <= cost_std && delta_count == 1
        lb = U_new_grad_points(i2);
        i_min = i2;
        delta_count = delta_count + 1;
    elseif abs(delta_obj(i2)) >= cost_std && delta_count == 2
        ub = U_new_grad_points(i2);
        delta_count = delta_count + 1;
        i_max = i2;
    end
end
% delta_count 
%in case no ub could be computed
if delta_count == 2
    ub = lb + 3;
    i_max = i_min + 4;
end

%determine maximum doe criterium values
if DoE_crit(i_min) >= DoE_crit(i_max)
    doe_out = lb-U_curr;
else
    doe_out = ub-U_curr;
end

lb_doe = lb;
ub_doe = ub;

end
