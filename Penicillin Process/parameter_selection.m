function [para_selection] = parameter_selection(U_vec,K_model,num_samples,...
    rel_range,num_selected_paras,case_study,sampling,num_outputs,sensitivity_crit,...
    initial_cond,i_b2b,extra_var,grad_dir_vec,grad_step_size_vec)
%function to perform a parameter selection with respect to model fitting
%(sensitivities of model ouput) as well as gradient correction
%(sensitivities of cost function and constraint gradients)

disp('PARAMETER SELECTION');

%number of parameters 
num_para = length(K_model);
        
%number of timepoints
num_sample_points = length(sampling);
         
%define range of deviation in the parameter values
parameter_range = rel_range*K_model;
        
%generate a latin hypercube design to sample the through the uncertain parameter space
X_sample = lhsdesign(num_samples,num_para);

%define model
model = @penicillin_process_model;

%define initial conditions
Y0 = initial_cond;
Y0_grad1 = Y0;
Y0_grad2 = Y0;

%define current operating point
U_curr = U_vec(i_b2b,:);

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

%initialize vectors and matrices      
para_selection = zeros(size(K_model));
Pe_matrix = [];
tc = [];
grad_sens_matrix = [];
% sens_local_sample_sum = zeros(num_para,1);
% sens_gradient_lhs = zeros(num_para,1);
% sens_gradient_lhs2 = zeros(num_para,1);
% sens_objective_sum = zeros(num_para,1);
% rank_non_full = 0;
% Pe_sample = [];
% Pe_sample2 = [];
% Pe_sample3 = [];
        
%for each sample of parameter values, calculate sensitivities   
for i_sample = 1:num_samples
    
    %get vector of samples
    X_sample_row = X_sample(i_sample,:);
    %from sample vector, get parameter values
    K_sample = getParaSample(K_model,parameter_range,X_sample_row);
    
    
    if ismember('output',sensitivity_crit)
        %obtain local sensitivities
        if strcmp(case_study,'penicillin')
            Y0(3) = U_curr(1);
            [~,output,local_sens] = sens_sys('penicillin_sens_analysis_model',sampling,Y0,opt,K_sample',[],[],U_curr,extra_var);
        end
        %NTxNYxNP

        i1=1;
        i2=1;
        i3=1;
        sens_local_y_p = zeros(num_outputs,num_para);
        sens_local_p = zeros(num_para,1);  
        output_magn = zeros(num_outputs,1);

        %calculate average output magnitude for each output
        output_magn = mean(abs(output));

            %get parameter estimiability matrix  
            for i_pe2 = 1:num_outputs
                for i_pe = 1:num_sample_points
                    %sensitivity with average output scaling
                    Pe_element(1:num_para) = local_sens(i_pe,i_pe2,:)./output_magn(i_pe2);

                    %parameter magnitude scaling 
                    %Already included in sensitivity
                    Pe_element = Pe_element.*K_sample;

                    %into matrix form
                    Pe_matrix = [Pe_matrix; Pe_element];
                end
            end
    end
    
    if ismember('gradient',sensitivity_crit)
        
        %calculate gradient sensitivities for each sample of the parameter
        %values
        
        
        %calculate nominal gradients
            if strcmp(case_study,'penicillin');
                Y0_grad1(3) = U_curr(1);
                [~, Y_nom] = ode15s(model,sampling, Y0_grad1, opt,U_curr,K_sample,extra_var);
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
                    [~, Y_grad] = ode15s(model,sampling, Y0_grad2, opt,U_grad,K_sample,extra_var);

                    %compute gradient
                    grad_pred_nominal(i_grad2) = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/...
                        (sqrt(grad_dev(1)^2+grad_dev(2)^2));  
                end
            end
            
        %calculate gradients with a deviation in parameter values
        for i_para = 1:num_para
            
            %define new parameter values
            K_new = K_sample;
            
            %deviation in parameter for sensitivity step size
            dev_para = sqrt(eps)*(K_new(i_para)+0.1);
            
            %change para value by small deviation
            K_new(i_para) = K_new(i_para)+dev_para;
            
            %calculate for nominal conditions
            if strcmp(case_study,'penicillin');
                Y0_grad1(3) = U_curr(1);
                [~, Y_nom] = ode15s(model,sampling, Y0_grad1, opt,U_curr,K_new,extra_var);
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
                    [~, Y_grad] = ode15s(model,sampling, Y0_grad2, opt,U_grad,K_new,extra_var);

                    %compute gradient
                    grad_pred_new(i_grad2) = (Y_grad(end,2)*Y_grad(end,4)-Y_nom(end,2)*Y_nom(end,4))/...
                        (sqrt(grad_dev(1)^2+grad_dev(2)^2));  
                end
            end
            
            %calculate gradient sensitivities
            grad_sens = (grad_pred_new-grad_pred_nominal)./abs(dev_para);
            
            %scale with nominal parameter and gradient value
            grad_sens =  grad_sens./abs(grad_pred_nominal);
            
            %put gradient sensitivities into matrix
            grad_sens_matrix(:,i_para) = grad_sens;
        end
        
    end

end

%initialize selection procedure
New_pe_matrix = Pe_matrix;
New_grad_sens_matrix = grad_sens_matrix;
for i_ortho = 1:min(num_selected_paras,length(K_model))
    
    %update tc
    tc_in = tc;
    
    %scale output sensitivities by their average magnitude
    if ismember('output',sensitivity_crit)
        %get magnitude of output sensitivities
        output_sens_ranking_vector = sum(abs(New_pe_matrix),1);

        %scale by average magnitude
        output_sens_ranking_vector = output_sens_ranking_vector./mean(output_sens_ranking_vector);
    else
        output_sens_ranking_vector = zeros(size(K_model));
    end

    %scale gradient sensitivities by their average magnitude
    if ismember('gradient',sensitivity_crit)
        %get magnitude of output sensitivities
        gradient_sens_ranking_vector = sum(abs(New_grad_sens_matrix));

        %scale by average magnitude
        gradient_sens_ranking_vector = gradient_sens_ranking_vector./mean(gradient_sens_ranking_vector);
    else
        gradient_sens_ranking_vector = zeros(size(K_model));
    end
    
    %overall sensitivity ranking vector
    overall_sens_ranking_vec = output_sens_ranking_vector + gradient_sens_ranking_vector;
    
    %determine parameter corresponding to maximum overall sensitivity
    [max_val,max_ind] = max(overall_sens_ranking_vec);
    
    %update parameter selection vector
    para_selection(max_ind) = 1;
    
    %perform orthogonalization to account for correlation between
    %parameters's effects on the model output
    if ismember('output',sensitivity_crit)
    [New_pe_matrix,tc] = perform_orthogonalization(Pe_matrix,max_ind,tc_in);
    end
    
    %update gradient sensitivity matrix
    New_grad_sens_matrix(:,max_ind) = zeros(size(grad_sens_matrix,1),1);
end




    function [K_sample] = getParaSample(K,range,X_row)

        %go through sample row
        for i_getpara1 = 1:length(K)
             %transform row value into paravalue
             K_sample(i_getpara1) = X_row(i_getpara1)*2*range(i_getpara1)+K(i_getpara1)-range(i_getpara1);
        end    
    end

    function [Res_pe_matrix,tc] = perform_orthogonalization(Original_pe_matrix,max_sens_para,tc)
        
        %parameter estimation matrix
        Pem = Original_pe_matrix;
         
        %get column of top_para
        tc = [tc Pem(:,max_sens_para)];
                
        %get lr matrix
        lr_matrix = tc*inv(tc'*tc)*tc'*Pem;
                
        %residual matrix
        res_matrix = Pem-lr_matrix;
        
        %new matrix used for ranking
        Res_pe_matrix = res_matrix;
        
    end
end