function [K_model_prime,new_corr,pred_grad] = gradient_correction(K_model,plant_gradient_vec,U_vec,i_b2b,...
    case_study,para_limit,lb,ub,initial_cond,para_selection,prev_corr,grad_dir_vec,...
    grad_step_size_vec,sampling,extra_var,num_corrected_grad,t_error,...
    batch_output_max_history,batch_output_min_history,num_outputs)
%function to perform gradient correction by adapting the model parameters

disp('GRADIENT CORRECTION STEP');

if strcmp(case_study,'penicillin')
    %Define penicillin process model
    model = @penicillin_process_model;
end

%Initial parameter values
K_ini = K_model;

%in case of further limitation in parameter values
if para_limit == 0
    lb_new = lb;
    ub_new = ub;
else
    rel_var = para_limit;
    lb_new=max(K_ini.*(1-rel_var),lb);
    ub_new=min(K_ini.*(1+rel_var),ub);
end

%current operating point
U_curr = U_vec(i_b2b,:);

%initial process conditions
Y0 = initial_cond;
Y0(3) = U_curr(1); 
Y0_grad1 = Y0;
Y0_grad2 = Y0;
Y0_old = Y0;



opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

%compute first order derivatives necessary for the correction term
[~,DX,DXDP] = sens_sys('penicillin_sens_analysis_model',sampling,Y0,opt,K_ini',[],[],U_curr,extra_var);


%fmincon options
options=optimset('Algorithm','active-set','Display','iter','UseParallel','always',...
    'TolX',1e-4,'TolFun',1e-4,'TolCon',1e-4,'MaxFunEvals',1500);

%'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-8

% if i_b2b < 2
% [x,fval] = fmincon(@gradobj,K_ini,[],[],[],[],lb_new,ub_new,[],options);
% else
    [x,fval] = fmincon(@gradobj,K_ini,[],[],[],[],lb_new,ub_new,@grad_constr,options);
% end

% @grad_constr
% @set_based_grad_constr1
% @set_based_grad_constr_prev_batch

K_model_prime = K_ini;
for i2 = 1:length(K_ini)
  if para_selection(i2) == 1
     K_model_prime(i2) = x(i2);
  end
end

%calculate correction term
Y_old_para = DX;
[~,Y_new_para] = ode15s(model,sampling, Y0, opt,U_curr,K_model_prime,extra_var);
c_new = -Y_old_para+Y_new_para;
new_corr=prev_corr+c_new;

%corrected gradients
pred_grad = grad_pred_output;



    function F2 = gradobj(x)
        
        K_new = K_ini;
        %only change the parameters which have been selected for gradient
        %correction
        
        for i_selection = 1:length(K_model)
          if para_selection(i_selection) == 1
              K_new(i_selection) = x(i_selection);
          end
        end
        
        %obtain errors between measured and predicted gradients
        %using the last "num_corrected_grad" gradient measurements 
        for i_corr_grad = 1:min(i_b2b,num_corrected_grad)
            
            %get correspodning operating point and 
            U_new = U_vec(i_b2b+1-i_corr_grad,:);
            
            plant_grad_new = plant_gradient_vec(i_b2b+1-i_corr_grad,:);
            
            %call up gradient error function
            F1(i_corr_grad) = compute_gradient_error(K_new,plant_grad_new,U_new);
           
        end
        
        %compute sum of gradient errors
        F2 = sum(F1.^(2));
%         F2 = sum(F1.^(2))+1e+2*norm(x-K_ini)^2;
        
     
    

    function F = compute_gradient_error(x,plant_grad,U)
         
        %calculate for nominal conditions
        Y0_grad1(3) = U(1);
        [~, Y_nom] = ode15s(model,sampling, Y0_grad1, opt,U,x,extra_var);
        

        %calculate predicted gradient
        for i_grad2 = 1:size(grad_dir_vec,1)

            %in case of penicillin process
            if strcmp(case_study,'penicillin')
                %update initial substrate concentration


                %update decision variables
                grad_dev = grad_step_size_vec(i_grad2)*grad_dir_vec(i_grad2,:);
                U_grad = U+grad_dev;

                Y0_grad2(3) = U_grad(1);
                [T_grad, Y_grad] = ode45(model,sampling, Y0_grad2, opt,U_grad,x,extra_var);

                %compute gradient
                grad_pred(i_grad2) = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/...
                    (sqrt(grad_dev(1)^2+grad_dev(2)^2));  
                
                if i_corr_grad == 1
                    grad_pred_output = grad_pred;
                end
            end
        end
        
        %obtain error between predicted and measured gradient
        
        %weighting matrix
%         W_grad = diag(plant_grad);    
        W_grad = diag([1 100]);  
%         W_grad = diag([1 10]);  
        
        F = sum(abs((plant_grad-grad_pred)/W_grad));
        
    end
    
    end

    function [myc,myceq] = grad_constr(x)
        
        K_constr = K_ini;
        %only change the parameters which have been selected for gradient
        %correction
        
        for i_selection = 1:length(K_ini)
          if para_selection(i_selection) == 1
              K_constr(i_selection) = x(i_selection);
          end
        end
    
       
       %calculate new model outputs
       if strcmp(case_study,'penicillin') 
        [~, Y_constr] = ode15s(model,sampling, Y0, opt,U_curr,K_constr,extra_var);
       end


        c = 0;
        for i_match = 1:length(K_ini)
            if para_selection(i_match) == 1
                c = c + DXDP(:,:,i_match)*(x(i_match)-K_ini(i_match));
            end
        end
        
% %                 c = 0;
% %         i_match_count = 1;
% %         for i_match = 1:length(K_ini)
% %             if para_selection(i_match) == 1
% %                 c = c + DXDP(:,:,i_match)*(x(i_match_count)-K_ini(i_match_count));
% %                 i_match_count = i_match_count + 1;
% %             end
% %         end

        Y_constr=Y_constr-c-prev_corr;  %new model
        Y_prime=DX-prev_corr;  %previous model
        trunc_e=abs(Y_constr./Y_prime-1);
        % trunc_e=abs(Y_prime./Y-1);



        % %avoid problem that parameter change due to gradient correction is too
        % %large in the second iteration
        % if (i == 2)
        %     myc=max(max(trunc_e))-t_error;
        % else
        %     %myc=max(max(trunc_e))-t_error;
        %     myc=max(trunc_e(end,:))-t_error;
        % end
        % t_error2 = noise_lvl;
        % t_error2 = 0.05;
        t_error2 = 0.25;
        % t_error2 = 5;
        % t_error2 = 10;
            
%         myc(1)= max(max(trunc_e(:,:)))-t_error;
        myc(1)= max(trunc_e(end,:))-t_error;  
        %test: relaxation of truncation error
        myc(2)= max(max(trunc_e(:,:)))-t_error2;

        myceq=[];

    end

 

end

