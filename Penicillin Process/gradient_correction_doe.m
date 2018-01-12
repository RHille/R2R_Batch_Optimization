function [grad_point_vec,max_doe_vec] = gradient_correction_doe(K_model,initial_cond,U_vec,sampling,extra_var,...
    U_out_vec,i_b2b,num_grad_selec,para_selection,num_grad_points)

%define model
model = @penicillin_process_model;

Y0 = initial_cond;
Y0_grad1 = Y0;
Y0_grad2 = Y0;
U_curr = U_vec(i_b2b,:);

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

num_para = length(K_model);

%calculate gradient sensitivities from output measurements from past
%operating points

%calculate nominal gradients
Y0_grad1(3) = U_curr(1);
[~, Y_nom] = ode15s(model,sampling, Y0_grad1, opt,U_curr,K_model,extra_var);


for i_grad_select = 1:min(i_b2b+1,num_grad_selec)

        %define difference in U
        %get deviation between correspodning operating points 
        U_new = U_out_vec(end-(min(i_b2b+1,num_grad_selec))+i_grad_select-1,:);
        U_dev_new = U_new-U_curr;
 
        Y0_grad2(3) = U_new(1);
        [~, Y_grad] = ode15s(model,sampling, Y0_grad2, opt,U_new,K_model,extra_var);

        %compute gradient
        grad_pred_nominal(i_grad_select) = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/...
        (sqrt(U_dev_new(1)^2+U_dev_new(2)^2));  
          
end
            

            
        %calculate gradients with a deviation in parameter values
        i_2 = 1;
        for i_para = 1:num_para
            
         %only calculate sensitivity for parameters which are used for
         %gradient correction
          if para_selection(i_para) == 1
              
            %define new parameter values
            K_new = K_model;
            
            %deviation in parameter for sensitivity step size
            dev_para = sqrt(eps)*(K_new(i_para)+0.1);
            
            %change para value by small deviation
            K_new(i_para) = K_new(i_para)+dev_para;
            
            %calculate for nominal conditions
           
                Y0_grad1(3) = U_curr(1);
                [~, Y_nom] = ode15s(model,sampling, Y0_grad1, opt,U_curr,K_new,extra_var);
          
            
            for i_grad_select = 1:min(i_b2b+1,num_grad_selec)

                    %define difference in U
                    %get deviation between correspodning operating points 
                    U_new = U_out_vec(end-(min(i_b2b+1,num_grad_selec))+i_grad_select-1,:);
                    U_dev_new = U_new-U_curr;

                    Y0_grad2(3) = U_new(1);
                    [~, Y_grad] = ode15s(model,sampling, Y0_grad2, opt,U_new,K_new,extra_var);

                    %compute gradient
                    grad_pred_new(i_grad_select) = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/...
                    (sqrt(U_dev_new(1)^2+U_dev_new(2)^2));  

            end
     
            %calculate gradient sensitivities
            grad_sens = (grad_pred_new-grad_pred_nominal)./abs(dev_para);
            
            %scale with nominal parameter and gradient value
%             grad_sens =  grad_sens./abs(grad_pred_nominal).*K_new(i_para);

%             grad_sens = grad_sens.;
            
            
            %put gradient sensitivities into matrix
            grad_sens_matrix(:,i_2) = grad_sens;
            
            i_2 = i_2 + 1;
          end
      end
      
%%%%%%%%%%%%%%%%%%EXPERIMENTAL DESIGN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DoE_crit = [];

grad_sens_matrix
U_out_vec
%use the two extra gradient measurements to obtain initial FIM
FIM_initial = grad_sens_matrix(end-1:end,:)'*grad_sens_matrix(end-1:end,:);

%get initial grad point vector
grad_point_vec = zeros(size(grad_sens_matrix,1),1);
grad_point_vec(end) = 1;
grad_point_vec(end-1) = 1;

%calculate doe objective in sequential fashion for number of selected
%gradient measurement points
FIM_element = grad_sens_matrix(end-1:end,:);

max_doe_vec = 0;
if i_b2b > 1
    for i_grad_points = 1:min(i_b2b,num_grad_points)

        %go through all operating points which are considered for gradient correction
        for i_grad_crit = 1:size(grad_sens_matrix,1)-2

            FIM_element_new = [FIM_element; grad_sens_matrix(i_grad_crit,:)];
            FIM_new = FIM_element_new'*FIM_element_new;
            DoE_crit_2(i_grad_crit) = det(FIM_new+0*eye(size(FIM_new)));  

%             FIM_element_new = [FIM_element; grad_sens_matrix(i_grad_crit,:)];
%              DoE_crit_2(i_grad_crit) = det(FIM_element+grad_sens_matrix(i_grad_crit,:)'*grad_sens_matrix(i_grad_crit,:));  
                 
        end
    
        %set values that have been selected to zero
        DoE_crit_3 = DoE_crit_2.*(grad_point_vec(1:end-2) == 0)';
               
        
% %         [~,a] = size(DoE_crit_3);
% %         
% %         
% %         t = fliplr(linspace(0,a,a));
% %         lambda = 0.85;
% %         a = 1 - 0.05*t;
% %         b = max(0,a);
        
%         size(DoE_crit_3)
%         size(t)
%         lambda.^((t))
%         DoE_crit_3
%         DoE_crit_3 = DoE_crit_3.*lambda.^((t));
%         DoE_crit_3 = DoE_crit_3.*b;
        
        %plot criterion
        figure(35+i_grad_points)
        scatter(U_out_vec(end-size(grad_sens_matrix,1):end-3,1)-U_out_vec(end,1),DoE_crit_3,'fill')
        drawnow

        %determine points which provides the maximum doe objective value
        [max_doe,ind_max] = max(DoE_crit_3);
        
        %get value of doe criterion
        max_doe_vec(i_grad_points) = max_doe;

        %update gradient point output vector
        grad_point_vec(ind_max) = 1

        %update FIM element
        FIM_element = [FIM_element; grad_sens_matrix(ind_max,:)];
%         FIM_element = [FIM_element'*FIM_element+grad_sens_matrix(ind_max,:)'*grad_sens_matrix(ind_max,:)];

    end
else
    grad_point_vec = [1;1];
    
end

% % %calculate the criterion for the remaining operating points
% % for i_grad_crit = 1:size(grad_sens_matrix,1)-2
% %     
% %     DoE_crit(i_grad_crit) = det(FIM_initial + grad_sens_matrix(i_grad_crit,:)'*grad_sens_matrix(i_grad_crit,:));
% %     
% % %     DoE_crit(i_grad_crit) = det(grad_sens_matrix(i_grad_crit,:)'*grad_sens_matrix(i_grad_crit,:));
% %     %test mit anderer Berechnung der next FIM
% %     FIM_element = [grad_sens_matrix(end-1:end,:); grad_sens_matrix(i_grad_crit,:)];
% %     DoE_crit_2(i_grad_crit) = det(FIM_element'*FIM_element);
% % end
% % DoE_crit'
% % %plot criterion for differnt operating points
% % if i_b2b > 1
% % figure(32)
% % scatter(U_out_vec(end-size(grad_sens_matrix,1):end-3,1)-U_out_vec(end,1),DoE_crit,'fill')
% % hold on
% % scatter(U_out_vec(end-size(grad_sens_matrix,1):end-3,1)-U_out_vec(end,1),DoE_crit_2,'k','fill')
% % hold off
% % % scatter(U_out_vec(end-size(grad_sens_matrix,1):end-1,1)-U_out_vec(end,1),DoE_crit,'fill')
% % drawnow
% % end
% % 
% % %test
% % %get max doe point and compute DoE_crit again
% % [~,ind_max] = max(DoE_crit);
% % FIM_element = [grad_sens_matrix(end-1:end,:); grad_sens_matrix(ind_max,:)];
% % for i_grad_crit = 1:size(grad_sens_matrix,1)-2
% %     
% %     DoE_crit_new(i_grad_crit) = det(FIM_initial + grad_sens_matrix(ind_max,:)'*grad_sens_matrix(ind_max,:) + grad_sens_matrix(i_grad_crit,:)'*grad_sens_matrix(i_grad_crit,:));
% %     
% %     FIM_element = [FIM_element; grad_sens_matrix(i_grad_crit,:)];
% %     DoE_crit_new_2(i_grad_crit) = det(FIM_element'*FIM_element);
% % %     DoE_crit(i_grad_crit) = det(grad_sens_matrix(i_grad_crit,:)'*grad_sens_matrix(i_grad_crit,:));
% % end
% % if i_b2b > 1
% % figure(33)
% % scatter(U_out_vec(end-size(grad_sens_matrix,1):end-3,1)-U_out_vec(end,1),DoE_crit_new,'r','fill')
% % hold on
% % scatter(U_out_vec(end-size(grad_sens_matrix,1):end-3,1)-U_out_vec(end,1),DoE_crit_new_2,'k','fill')
% % hold off
% % % scatter(U_out_vec(end-size(grad_sens_matrix,1):end-1,1)-U_out_vec(end,1),DoE_crit,'fill')
% % drawnow
% % end

end