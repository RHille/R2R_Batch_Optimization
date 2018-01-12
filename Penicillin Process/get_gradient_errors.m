function [grad_error_vector] = get_gradient_errors(K_model,U_vec,i_b2b,...
    case_study,initial_cond,para_selection,sampling,extra_var,U_out_vec,Obj_fun_vec,grad_point_vec)

    %function to calculate the errors between measured and predicted
    %gradients and normalized cost function differences
    
    if strcmp(case_study,'penicillin')
        %Define penicillin process model
        model = @penicillin_process_model;
    end
    
    %current operating point
    U_curr = U_vec(i_b2b,:);

    %initial process conditions
    Y0 = initial_cond;
    Y0(3) = U_curr(1); 
    Y0_grad1 = Y0;
    Y0_grad2 = Y0;
    Y0_old = Y0;
    
    
    %Initial parameter values
    K_ini = K_model;
    
    
    opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
    
    
    F2 = gradobj(K_model);
    
    grad_error_vector = flipud(F1');
    
    
    function F2 = gradobj(x)
        
        plant_grad_vec = [];
        grad_pred_vec = [];
        
        K_new = K_ini;
        %only change the parameters which have been selected for gradient
        %correction
        
        for i_selection = 1:length(K_model)
          if para_selection(i_selection) == 1
              K_new(i_selection) = x(i_selection);
          end
        end
        
        
        %calculate for nominal conditions
        U_curr = U_out_vec(end,:);
        Y0_grad1(3) = U_curr(1);
        [~, Y_nom] = ode15s(model,sampling, Y0_grad1, opt,U_curr,K_new,extra_var);
        
        %obtain errors between measured and predicted gradients
        %using the last "num_corrected_grad" gradient measurements 
        for i_corr_grad = 1:length(grad_point_vec)
            
            %compute gradient if operating point has been selected for
            %gradient correction
%             if grad_point_vec(end-i_corr_grad+1) == 1
                %get deviation between correspodning operating point and 
                U_new = U_out_vec(end-i_corr_grad,:);
                U_dev_new = U_new-U_curr;

                %calculate difference between cost measurements
                Cost_dev_new = Obj_fun_vec(end-i_corr_grad)-Obj_fun_vec(end);

                plant_grad_new = Cost_dev_new/sqrt(sum(U_dev_new.^2));

                %call up gradient error function
                F1(i_corr_grad) = compute_gradient_error(K_new,plant_grad_new,U_curr,U_new,U_dev_new);
%             end
           
        end
        
        %compute sum of gradient errors
        F2 = sum(F1.^2);
        
     
    

    function F = compute_gradient_error(x,plant_grad,U_curr,U_new,U_dev_new)
        
        %calculate predicted gradient
%         for i_grad2 = 1:size(grad_dir_vec,1)

            %in case of penicillin process
            if strcmp(case_study,'penicillin')
                %update initial substrate concentration


                %update decision variables
%                 grad_dev = grad_step_size_vec(i_grad2)*grad_dir_vec(i_grad2,:);
%                 U_grad = U+grad_dev;

                Y0_grad2(3) = U_new(1);
                [T_grad, Y_grad] = ode15s(model,sampling, Y0_grad2, opt,U_new,x,extra_var);

                %compute gradient
                grad_pred = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/sqrt(sum(U_dev_new.^2));  
                
%                 if i_corr_grad == 1
%                     grad_pred_output = grad_pred;
%                 end
            end
%         end
        
        plant_grad_vec = [plant_grad_vec plant_grad];
        grad_pred_vec = [grad_pred_vec grad_pred];
        
        %obtain error between predicted and measured gradient
        
        %weighting matrix
        W_grad = diag(plant_grad);    
%         W_grad = diag([1 100]);  
        
        F = (abs((plant_grad-grad_pred)/W_grad));
%         F = (abs((plant_grad-grad_pred)));
        
    end
    
    end
    

end
