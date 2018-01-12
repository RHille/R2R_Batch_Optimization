function plot_grad_sens_objective_function(K_model,sampling,initial_cond,U_vec,para_selection,...
    i_b2b,proc_out_mean_history,num_outputs,extra_var,prev_corr,case_study,grad_step_size_vec,...
    grad_dir_vec)



%determine range for parameters

K_model_new = K_model;
K_ini = K_model_new;
count_limit = 10;

U_curr = U_vec(i_b2b,:);
Y0 = initial_cond;
%100% range
Y0_grad1 = Y0;
Y0_grad2 = Y0;

model = @penicillin_process_model;

%ode options
opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

proc_out_mean = reshape(proc_out_mean_history(:,i_b2b),length(sampling),num_outputs);

i_count3 = 1;
        for i = 1:length(para_selection)
            if para_selection(i) == 1
                K_model_vec(i_count3,:) = 0.01*K_model(i):0.1*K_model(i):2.01*K_model(i);
                i_count3 = i_count3 + 1;
            end
        end

%test for different operating points
% U_vec_new = [50:5:55];
U_vec_new = [55];
% K8_vec = [0.014:0.003:0.016];
U_vec_new = [0.015:0.005:0.018];
% for i_K8 = 1:length(K8_vec)
for i_u = 1:length(U_vec_new)
    
%     U_curr(1) = U_vec_new(i_u);
    U_curr(1) = 54;

K2_vec = 0.01:0.05:0.7;
K5_vec = [K2_vec*0.1];

K1_vec = 0.08:0.005:0.12;
K3_vec = 0.001:0.005:0.08;
i1 = 1;

for i_count = 1:length(K1_vec)
%     length(K2_vec)
    i2 = 1;
    for i_count2 = 1:length(K3_vec)
%         length(K5_vec)
        
        selec_count = 0;
%         for i = 1:length(para_selection)
%             if para_selection(i) == 1 && selec_count == 0
%                 K_model_new(i) = 0.01*K_model(i)+i_count/count_limit*K_model(i);
%                 selec_count = selec_count + 1;
%             elseif para_selection(i) == 1 && selec_count == 1
%                 K_model_new(i) = 0.01*K_model(i)+i_count2/count_limit*K_model(i);
%             end
% %         end
%         K_model_new(2) = K2_vec(i_count);
%         K_model_new(5) = K5_vec(i_count2);
        
        K_model_new(1) = K1_vec(i_count);
        K_model_new(3) = K3_vec(i_count2);
        K_model_new(8) = U_vec_new(i_u);
        
        
        F2 = grad_sens_obj(K_model_new); 
        
        grad_sens_cost_mat(i2,i1) = F2;
        
        i2 = i2 + 1;
    end
    i1 = i1 +1;
    disp('.');
end
if i_u == 1
    grad_sens_cost_mat_new = zeros(size(grad_sens_cost_mat));
end

grad_sens_cost_mat_new = grad_sens_cost_mat_new + min(grad_sens_cost_mat,45)/length(U_vec_new);
% grad_sens_cost_mat_new = grad_sens_cost_mat_new + grad_sens_cost_mat/length(U_vec_new);
disp('......');
end


figure(30)
    %contour plot
% %     size(K_model_vec(2,:))
% %     size(grad_sens_cost_mat)
    mean(mean(grad_sens_cost_mat))
    
    grad_sens_cost_mat2 = (grad_sens_cost_mat>7);
    grad_sens_cost_mat = min(grad_sens_cost_mat,mean(mean(grad_sens_cost_mat)));
%     contourf(K2_vec,K5_vec,grad_sens_cost_mat2)
%     surf(K2_vec,K5_vec,grad_sens_cost_mat)
%     surf(K2_vec,K5_vec,grad_sens_cost_mat_new)
%     contourf(K2_vec,K5_vec,grad_sens_cost_mat_new)
    contourf(K1_vec,K3_vec,grad_sens_cost_mat_new)
%     surf(K1_vec,K3_vec,grad_sens_cost_mat)
    xlabel('Gln');
    ylabel('Glc');
    zlabel('SSE');
    drawnow
    
    function F = grad_sens_obj(x)
            
            i_count5 = 1;
            K_new = K_ini;
             for i5 = 1:length(K_ini)
                if para_selection(i5) == 1
                    K_new(i5) = x(i5);
%                     K_new(i5) = x(i_count5);
%                     i_count5 = i_count5 + 1;
                end
             end
            
            %%calculate nominal gradients
            if strcmp(case_study,'penicillin');
                Y0_grad1(3) = U_curr(1);
                [~, Y_nom] = ode45(model,sampling, Y0_grad1, opt,U_curr,K_new,extra_var);
            end

            %calculate predicted gradient
            for i_grad2 = 1:size(grad_dir_vec,1)

                %in case of penicillin process
                if strcmp(case_study,'penicillin')
                    %update initial substrate concentration


                    %update decision variables
                    grad_dev = grad_step_size_vec(i_grad2)*grad_dir_vec(i_grad2,:);
                    U_grad = U_curr+grad_dev;

                    Y0_grad2(3) = U_grad(1);
                    [~, Y_grad] = ode45(model,sampling, Y0_grad2, opt,U_grad,K_new,extra_var);

                    %compute gradient
                    grad_pred_nominal(i_grad2) = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/...
                        (sqrt(grad_dev(1)^2+grad_dev(2)^2));  
                end
            end
        
        i6 = 1;
         %calculate gradients with a deviation in parameter values
        for i_para = 1:length(K_ini)
            if para_selection(i_para) == 1
            
            %define new parameter values
            K_new2 = K_new;
            
            %deviation in parameter for sensitivity step size
            dev_para = sqrt(eps)*(K_new2(i_para)+0.1);
            
            %change para value by small deviation
            K_new2(i_para) = K_new2(i_para)+dev_para;
            
            %calculate for nominal conditions
            if strcmp(case_study,'penicillin');
                Y0_grad1(3) = U_curr(1);
                [~, Y_nom] = ode45(model,sampling, Y0_grad1, opt,U_curr,K_new2,extra_var);
            end

            %calculate predicted gradient
            for i_grad2 = 1:size(grad_dir_vec,1)

                %in case of penicillin process
                if strcmp(case_study,'penicillin')
                    %update initial substrate concentration


                    %update decision variables
                    grad_dev = grad_step_size_vec(i_grad2)*grad_dir_vec(i_grad2,:);
                    U_grad = U_curr+grad_dev;

                    Y0_grad2(3) = U_grad(1);
                    [~, Y_grad] = ode45(model,sampling, Y0_grad2, opt,U_grad,K_new2,extra_var);

                    %compute gradient
                    grad_pred_new(i_grad2) = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/...
                        (sqrt(grad_dev(1)^2+grad_dev(2)^2));  
                end
            end
            
            %calculate gradient sensitivities
            grad_sens = abs(grad_pred_new-grad_pred_nominal)./abs(dev_para);
            
            %scale with nominal parameter and gradient value
            grad_sens =  grad_sens./abs(grad_pred_nominal).*K_new2(i_para);
%             grad_sens =  grad_sens./abs(grad_pred_nominal);
            
            %put gradient sensitivities into matrix
            grad_sens_matrix(:,i6) = grad_sens;
            
            i6 = i6 + 1;
            end
        end
        %test
%         grad_sens_matrix = grad_sens_matrix(:,1);
       
        F = sum(sum(grad_sens_matrix));

end


end