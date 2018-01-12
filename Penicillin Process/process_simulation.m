function [outputs_mean,outputs_max,outputs_min,grad_additive,grad_noisy,U_out_vec,Obj_fun_vec,pen_std] = process_simulation(sampling_time,number_batches,...
    measurement_noise_lvl,K_process,para_noise,initial_cond_uncertainty,...
    initial_cond,decision_var,extra_var,case_study,grad_step_size_vec,...
    grad_dir_vec,num_grad_batches,grad_noise_lvl,i_b2b,smooth_set_constr,...
    U_out_vec,Obj_fun_vec)
%function to generate process measurements used in the B2B optimization
%procedure

disp('PROCESS SIMULATION');

pen_simulator = @penicillin_process_simulator;

%current operating point
U_nom = decision_var(i_b2b,:);

sampling = sampling_time;
opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

initial_cond(3) = U_nom (1);
U = U_nom;

num_outputs = length(initial_cond);
%for specified number of batches, simulate process
for i_batch_run = 1:number_batches
    
    %Implement possible uncertainties in initial conditions
    
    if i_b2b > 1
        Y0=mvnrnd(initial_cond,diag(initial_cond_uncertainty));
    else
        Y0 = initial_cond;
    end
    
    %Implement possible uncertainties in process parameters (to simulate
    %process noise)
    K_sim = K_process.*(1+para_noise*randn(size(K_process)));
    
    %in case of penicillin process
    if strcmp(case_study,'penicillin')
        %update initial substrate concentration
        
        
%         Y0(3) = U(1);
        U(1) = Y0(3);
        [T, Y] = ode15s(pen_simulator,sampling, Y0, opt,U,K_sim,extra_var);
    end
    
    %create additive noise based on the average output magnitude 
    if i_batch_run == 1
        
        %calculate average magnitude
        y_avrg = mean(Y);
        
        %use 10% of the magnitude x the noise lvl as the standard deviation
        %for noise_lvl of 1 we get 10% measurement noise
        std_data = measurement_noise_lvl*0.1*y_avrg;
        
        %for set based constraints
        constr_data = std_data;
        if strcmp(case_study,'penicillin')
            %increase noise for penicillin measurement
            %test: relaxation of set-based constraints
             constr_data(2) = 3*constr_data(2);

        
        end
        %adjustment for penicillin
%         std_data(2) = measurement_noise_lvl*0.02;
%          std_data(2) = std_data(2);
         pen_std = 1*std_data(2)*y_avrg(4);
    end
    
    exp_data = [];
    exp_data2 = [];
    Y_run_upper = [];
    Y_run_lower = [];
    %generate measurement noise
    for i_noise = 1:num_outputs
        new_element = normrnd(Y(:,i_noise),std_data(i_noise),size(Y,1),1);
        exp_data = [exp_data new_element];
        
        if i_noise == 2
            exp_data2 = [exp_data2 normrnd(Y(:,i_noise),constr_data(i_noise),size(Y,1),1)];
        else
            exp_data2 = [exp_data2 new_element];
        end
        
        if strcmp(smooth_set_constr,'on');
            Y_run_upper = [Y_run_upper Y(:,i_noise)+1.96*ones(length(T),1)*constr_data(i_noise)];
            Y_run_lower = [Y_run_lower Y(:,i_noise)-1.96*ones(length(T),1)*constr_data(i_noise)];
% %             Y_run_upper = [Y_run_upper Y(:,i_noise)+1.26*ones(length(T),1)*constr_data(i_noise)];
% %             Y_run_lower = [Y_run_lower Y(:,i_noise)-1.26*ones(length(T),1)*constr_data(i_noise)];
        end
    end
    
    %transform into vector
    exp_data_mat(:,i_batch_run) = exp_data(:).*(exp_data(:)>0);
    exp_data_mat2(:,i_batch_run) = exp_data2(:);
    Y_data_mat(:,i_batch_run) = Y(:); 
    if strcmp(smooth_set_constr,'on');
        Y_run_upper_mat(:,i_batch_run) = max(Y_run_upper(:),-0.5*ones(length(Y_run_upper(:)),1));
        Y_run_lower_mat(:,i_batch_run) = max(Y_run_lower(:),-0.5*ones(length(Y_run_upper(:)),1));
    end
end


%get range of values for each output (max,min and mean)
for i_out_val = 1:num_outputs
    output_mat = exp_data_mat(1+(i_out_val-1)*length(T):i_out_val*length(T),:);
    outputs_mean(:,i_out_val) = mean(output_mat,2);
    
    if strcmp(smooth_set_constr,'on');
        upper_mat = Y_run_upper_mat(1+(i_out_val-1)*length(T):i_out_val*length(T),:);
        lower_mat = Y_run_lower_mat(1+(i_out_val-1)*length(T):i_out_val*length(T),:);
        outputs_max(:,i_out_val) = max(upper_mat,[],2);
        outputs_min(:,i_out_val) = min(lower_mat,[],2);
    else
        outputs_max(:,i_out_val) = max(output_mat,[],2);
        outputs_min(:,i_out_val) = min(output_mat,[],2);
    end
    
    %for noise free outputs
    Y_mat = Y_data_mat(1+(i_out_val-1)*length(T):i_out_val*length(T),:);
    Y_max(:,i_out_val) = max(Y_mat,[],2);
    Y_min(:,i_out_val) = min(Y_mat,[],2);
    Y_mean(:,i_out_val) = mean(Y_mat,2);

end



%Perform gradient measurements

    
%go through all gradient directions
for i_grad2 = 1:size(grad_dir_vec,1)
        
   for i_grad = 1:num_grad_batches
       
        %in case of penicillin process
        if strcmp(case_study,'penicillin')
            %update initial substrate concentration
            
            %Implement possible uncertainties in initial conditions
            if i_b2b > 1
                Y0=mvnrnd(initial_cond,diag(initial_cond_uncertainty));
            else
                Y0 = initial_cond;
            end
    
            %Implement possible uncertainties in process parameters (to simulate
            %process disturbances)
            K_sim = K_process.*(1+para_noise*randn(size(K_process)));
            
            %update decision variables
            grad_dev = grad_step_size_vec(i_grad2)*grad_dir_vec(i_grad2,:);
            U(1) = Y0(3);
            U_grad = U+grad_dev;
            
            Y0(3) = U_grad(1);
            [T, Y] = ode15s(pen_simulator,sampling, Y0, opt,U_grad,K_sim,extra_var);
            
            
            exp_data_grad = [];
            %generate measurement noise
            for i_noise = 1:num_outputs
                exp_data_grad = [exp_data_grad normrnd(Y(:,i_noise),std_data(i_noise),size(Y,1),1)];
            end
    
            %transform into vector
            exp_data_grad_mat(:,i_grad) = exp_data_grad(:);
            Y_grad_mat(:,i_grad) = Y(:);
        end
        
       %get average values from gradient measurements
       for i_out_val = 1:num_outputs
           
           output_grad_mat = exp_data_grad_mat(1+(i_out_val-1)*length(T):i_out_val*length(T),:);
           output_grad_max(:,i_out_val) = max(output_grad_mat,[],2);
           output_grad_min(:,i_out_val) = min(output_grad_mat,[],2);
           output_grad_mean(:,i_out_val) = mean(output_grad_mat,2);
           
           Y_grad = Y_grad_mat(1+(i_out_val-1)*length(T):i_out_val*length(T),:);
           Y_grad_max(:,i_out_val) = max(Y_grad,[],2);
           Y_grad_min(:,i_out_val) = min(Y_grad,[],2);
           Y_grad_mean(:,i_out_val) = mean(Y_grad,2);
       end   
   end
   
       if strcmp(case_study,'penicillin')
%         if i_grad2 == 1
       %Save objective function values at corresponding points of the decision
        %variable
            U_out_vec = [U_out_vec; U_nom+grad_dev];
            Obj_fun_vec = [Obj_fun_vec; output_grad_mean(end,2)*Y_grad_mean(end,4)];
%         end
       end
  
   %calculate mean gradient
   if strcmp(case_study,'penicillin')
       
       %calculate noise free gradient
       grad_noise_free(i_grad2) = (Y_grad_mean(end,2)*Y_grad_mean(end,4)-Y_mean(end,2)*Y_mean(end,4))...
           /(sqrt(grad_dev(1)^2+grad_dev(2)^2));   
       
       %calculate gradient from noisy measurements
       grad_noisy_measurements(i_grad2) = (output_grad_mean(end,2)*Y_grad_mean(end,4)-outputs_mean(end,2)*Y_mean(end,4))...
           /(sqrt(grad_dev(1)^2+grad_dev(2)^2));   
   end
   
end

%Save objective function values at corresponding points of the decision
%variable
if strcmp(case_study,'penicillin') 
    U_out_vec = [U_out_vec; U_nom];
    Obj_fun_vec = [Obj_fun_vec; outputs_mean(end,2)*Y_mean(end,4)];
end

if strcmp(case_study,'penicillin')
    
    %estimated offline
    avrg_grad_mag = [2.664;123];
%     avrg_grad_mag = [2.664;123;2.664;2.664];
    
    %put additive gaussian noise on noise free gradient measurements
%     grad_additive = mvnrnd(grad_noise_free',grad_noise_lvl*diag(avrg_grad_mag));
    grad_additive = grad_noise_free;
    
    %noisy gradients
    grad_noisy = grad_noisy_measurements;
end

%relaxation of set-based constraints
[out_max,out_min] = get_sb_constraints(exp_data_mat2);

%update set based constraints
outputs_max = out_max;
outputs_min = out_min;


    function [out_max,out_min] = get_sb_constraints(exp_data_mat)
        
        %get max and min output values
        for i_out_val2 = 1:num_outputs
            output_mat = exp_data_mat(1+(i_out_val2-1)*length(T):i_out_val2*length(T),:);
            
            max_output(:,i_out_val2) = max(output_mat,[],2);
            min_output(:,i_out_val2) = min(output_mat,[],2);
            mean_output(:,i_out_val2) = mean(output_mat,2);
            
            %implement central moving average filter for constraint
            %smoothing
            [y_max,y_min] = central_moving_average(T,max_output(:,i_out_val2)',min_output(:,i_out_val2)');
            
            [y_max2,y_min2] = moving_average(T,max_output(:,i_out_val2)',min_output(:,i_out_val2)');
           
            y_max_cma(:,i_out_val2) = y_max;
            y_min_cma(:,i_out_val2) = y_min;
            
            y_max_ma(:,i_out_val2) = y_max2;
            y_min_ma(:,i_out_val2) = y_min2;
            
            %simple smoothing factor
            [y_max3,y_min3] = smoothing_factor(T,mean_output(:,i_out_val2)',max_output(:,i_out_val2)',min_output(:,i_out_val2)');
            
            y_max_sf(:,i_out_val2) = y_max3;
            y_min_sf(:,i_out_val2) = y_min3;
        end
        
        %prevent too restrictive constraints for biomass
        y_max_cma(:,1) = max(y_max_cma(:,1),1);

        
        %add smoothing based on standard deviation
        y_max_cma(:,1) = max(y_max_cma(:,1),mean_output(:,1)+constr_data(1)*1.5);
        y_max_cma(:,2) = max(y_max_cma(:,2),mean_output(:,2)+constr_data(2)*2.5);
        y_max_cma(:,3) = max(y_max_cma(:,3),mean_output(:,3)+constr_data(3));
        
        y_min_cma(:,1) = min(y_min_cma(:,1),mean_output(:,1)-constr_data(1)*1.5);
        y_min_cma(:,2) = min(y_min_cma(:,2),mean_output(:,2)-constr_data(2)*2.5);
        y_min_cma(:,3) = min(y_min_cma(:,3),mean_output(:,3)-constr_data(3));
        
        
%         out_max = y_max_cma;
%         out_min = y_min_cma;

        y_min_sf = max(y_min_sf,-0.1);
        
        out_max = y_max_sf;
        out_min = y_min_sf;
        
% %         %test for noise-free case
% %         out_max = y_max_sf+5;
% %         out_min = y_min_sf-5;
        
%         %plot max and min output values
%         figure(30)
%         clf
%         hold on
%         for i_out_val2 = 1:num_outputs-1
%             output_mat = exp_data_mat(1+(i_out_val2-1)*length(T):i_out_val2*length(T),:);
%             
%             for i2 = 1:size(output_mat,2)
%                 scatter(T,output_mat(:,i2),'b')
%             end
%         end
       
        
        
        
        
%         grid on;
%         plot(T,max_output(:,1),'b')
%         plot(T,min_output(:,1),'b')
% %         plot(T,y_max_cma(:,1),'--b')
% %         plot(T,y_min_cma(:,1),'--b')
% %         plot(T,outputs_max(:,1),':b')
% %         plot(T,outputs_min(:,1),':b')
% %         plot(T,y_max_ma(:,1),':b')
% %         plot(T,y_min_ma(:,1),':b')
%         plot(T,y_max_sf(:,1),'--b')
%         plot(T,y_min_sf(:,1),'--b')
%         plot(T,max_output(:,2),'r')
%         plot(T,min_output(:,2),'r')
% %         plot(T,y_max_cma(:,2),'--r')
% %         plot(T,y_min_cma(:,2),'--r')
% %         plot(T,outputs_max(:,2),':r')
% %         plot(T,outputs_min(:,2),':r')
% %         plot(T,y_max_ma(:,2),':r')
% %         plot(T,y_min_ma(:,2),':r')
%         plot(T,y_max_sf(:,2),'--r')
%         plot(T,y_min_sf(:,2),'--r')
%         plot(T,max_output(:,3),'k')
%         plot(T,min_output(:,3),'k')
% %         plot(T,y_max_cma(:,3),'--k')
% %         plot(T,y_min_cma(:,3),'--k')
% %         plot(T,outputs_max(:,3),':k')
% %         plot(T,outputs_min(:,3),':k')
% %         plot(T,y_max_ma(:,3),':k')
% %         plot(T,y_min_ma(:,3),':k')
%         plot(T,y_max_sf(:,3),'--k')
%         plot(T,y_min_sf(:,3),'--k')
%         drawnow
%         hold off
        
        function [y_max_cma,y_min_cma] = central_moving_average(sampling,max_output,min_output)
            
            y_max_cma = [];
            y_min_cma = [];
            
            y_max_cma(1) = max_output(:,1);
            y_min_cma(1) = min_output(:,1);
            
            for i=1:length(sampling)
            
                %previous
%                 y_max_nom_prev = max_output(:,i-min(1,i-1));
%                 y_min_nom_prev = min_output(:,i-min(1,i-1));
                y_max_nom_prev = y_max_cma(:,i-min(1,i-1));
                y_min_nom_prev = y_min_cma(:,i-min(1,i-1));


                %current
                y_max_nom_curr = max_output(:,i);
                y_min_nom_curr = min_output(:,i);

                %future
                y_max_nom_fut = max_output(:,i+min(1,length(sampling)-i));
                y_min_nom_fut = min_output(:,i+min(1,length(sampling)-i));

                %calculate central ma
                y_max_scma(i) = 1/3*(y_max_nom_prev+y_max_nom_curr+y_max_nom_fut);
                y_min_scma(i) = 1/3*(y_min_nom_prev+y_min_nom_curr+y_min_nom_fut);

                %use max min criteria
                y_max_cma(i) = max(y_max_scma(i),y_max_nom_curr);
                y_min_cma(i) = min(y_min_scma(i),y_min_nom_curr);

            end
            
            
        end
        
        function [y_max_ma,y_min_ma] = moving_average(sampling,max_output,min_output)
            
            y_max_ma = [];
            y_min_ma = [];
            y_max_sma = [];
            y_min_sma = [];

            %order
            n_ma = 2;
            lambda = 1;

            for i=1:length(sampling)

                %nominal max and min values
                y_max_nom(i) = max_output(:,i);
                y_min_nom(i) = min_output(:,i);
    %             y_max_ma(i) = y_max_nom(i);
    %             y_min_ma(i) = y_min_nom(i);

                %use ma filter
                max_sum = 0;
                min_sum = 0;
                for j=0:n_ma-1 
                    max_sum = max_sum+(lambda^min(i-1,j))*y_max_nom(i-min(i-1,j));
                    min_sum = min_sum+(lambda^min(i-1,j))*y_min_nom(i-min(i-1,j));

    %                 max_sum = max_sum+(lambda^min(i-1,j))*y_max_ma(i-min(i-1,j));
    %                 min_sum = min_sum+(lambda^min(i-1,j))*y_max_ma(i-min(i-1,j));
                end

                %update constraints
                y_max_sma(i) = 1/n_ma*max_sum;
                y_min_sma(i) = 1/n_ma*min_sum;

                %use max min criteria
                y_max_ma(i) = max(y_max_sma(i),y_max_nom(i));
                y_min_ma(i) = min(y_min_sma(i),y_min_nom(i));

            end
            
        end
        
        function [y_max_sf,y_min_sf] = smoothing_factor(sampling,mean_output,max_output,min_output)
        
        %calculate a simple smoothing factor to add to the mean
        y_max_factor(1) = 0;
        y_min_factor(1) = 0;
        y_mean(1) = mean_output(:,1);
        y_max_nom_curr(1) = max_output(:,1);
        y_min_nom_curr(1) = min_output(:,1);
        
        for i=2:length(sampling)
            
            %calculate mean
            y_mean(i) = mean_output(:,i);
            y_max_nom_curr(i) = max_output(:,i);
            y_min_nom_curr(i) = min_output(:,i);
            
            %difference of mean to max and min vals
            y_max_factor(i) = max(y_max_factor(i-1),y_max_nom_curr(i)-y_mean(i));
            y_min_factor(i) = max(y_min_factor(i-1),y_mean(i)-y_min_nom_curr(i));  
        end
        
        y_max_sf = y_mean+max(y_max_factor);
        y_min_sf = y_mean-max(y_min_factor);
        
        end
        
    end
end