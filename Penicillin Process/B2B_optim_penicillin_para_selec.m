function B2B_optim_penicillin_para_selec
%Main file to run the batch-to-batch optimization procedure of the
%penicillin process

%Thesis, Chapter3, Parameter Selection

%Author: Rubin Hille


%Program Initialization 
clc
close all
clear

%define case studies in list
case_study = {'penicillin'};
%% Penicillin Process Description


%Define penicillin process simulator
pen_simulator = @penicillin_process_simulator;

%process parameters
K_process=[0.092 0.15 0.005 0.0002 0.1 0.04 0.45 0.9 0.014];

%Define penicillin process simulator
pen_model = @penicillin_process_model;

%Initial process model parameters
K_model_initial=dlmread('fmincon_parameters_LCLS_2.txt');
%adjustment of parameters


%Bounds on model parameters (bounded parameter space)
lb_model_parameter = 6.5e-3*ones(size(K_model_initial));
ub_model_parameter = 5*ones(size(K_model_initial));


%extra bound for parameter 2
lb_model_parameter(2) = 2e-1;
ub_model_parameter(5) = 1e-1;
%extra bound for parameter 1
ub_model_parameter(1) = 0.1;

%Initial sampling time
sampling = [0:6:8*24];

%number of outputs
num_outputs = 4;

%list of output descriptions
list_output = {'biomass','penicillin','substrate','volume'};

%Initial substrate concentration
S0_initial = 0.1;

%Initial flowrate
F_initial = 0.01;

%Initial volume
vol_initial = 100;

%Initial biomass
bio_initial = 0.1;

%Initial penicillin
pen_initial = 0;

%Initial process conditions
Y_initial = [bio_initial pen_initial S0_initial vol_initial];

%Uncertainty in initial conditions in %
Y_initial_uncertainty = [0.0007 0 .5 0]*.8;

%In case of no input uncertainty
Y_initial_uncertainty = [0.0007 0 .5 0]*.0;

%Initial decision variables
U_initial = [S0_initial F_initial];

%Bounds on the decision variables
lb_opt = [0.1 0.01];
ub_opt = [100 100];

%Inlet Feed concentration
sf_inlet = 600;

%Constraint on the reactor volum
reactor_vol_constr = 120;

%uncertainty in process parameters to simulate process disturbances
para_uncetainty_lvl = 0;

%Define extra variables
extra_var = [sf_inlet];

%noise lvl in case that noise is added on noise free gradient measurements
%in [%]
grad_noise_lvl = 0.01;

P_final_process = [];
U_doe_vector = [];
%% Penicillin Process Optimization and Illustration

%For test purposes, find optimum and plot objective function for the "true"
%penicillin process

%call up function find optimum
% pen_process_optimum = optimize_penicillin_process(Y_initial,U_initial,sampling,extra_var,reactor_vol_constr,K_process);

%optimum of the actual process (process simulator w/o model-plant mismatch)
pen_process_optimum = [54.7338 0.1728];

%create a surface and contour plot of cost function of the actual process
% plot_penicillin_cost_function(pen_process_optimum,Y_initial,U_initial,extra_var,sampling,K_process)


%% B2B Algorithm Options

%upper bound on the relative truncation error

%number of batches per operating point
num_batch_per_op = 1;

%number of batches run for gradient measurements
num_grad_batches = 1;

%define level of measurement noise in the system. The level 1 corresponds
%to 10% measurement noise, level 0.5 would correspond to 5% additive
%gaussian noise
measurement_noise_lvl = 1;
% measurement_noise_lvl = 0.5;

%upper bound on the relative truncation error
eps_bound = 0.01;

%% B2B Algorithm Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define maximum number of batch runs
num_iterations = 35;

for i_sim = 1:10
    
%Information for gradient measurements
%gradient directions in terms of decision variables
grad_dir_vec = [1 0;0 1];


%step size in the gradient directions
grad_step_size_vec = [3; 0.5];


%Initial model parameters
K_model = K_model_initial;

%Initial operating point
U_new = [.1 0.04];
% U_new = [.1 0.1728];

%Vector of trajectory of decision variables
U_vec = [];
U_vec(1,:) = U_new;

%set initial correction term to zero
prev_corr = 0;
new_corr_mat = [];

%For storing of objective function measurements
U_out_vec = [];
Obj_fun_vec = [];

%for storing parameter values
para_estim_vec = [];

lb_doe = [];
ub_doe = [];

%start B2B procedure
for i_b2b = 1:num_iterations
%     
    
    %% %%%%%%%%%%%%%%%%% Step 1: Process simulation %%%%%%%%%%%%%%%%%%%%%%%

    %for smoothing of set-based constaints
    smooth_set_constr = 'on';
    
    %Call up function process simulation to simulate process
    [proc_out_mean,proc_out_max,proc_out_min,grad_additive,grad_noisy,...
        U_out_vec,Obj_fun_vec,pen_std] = process_simulation(sampling,...
        num_batch_per_op,measurement_noise_lvl,K_process,para_uncetainty_lvl,...
        Y_initial_uncertainty,Y_initial,U_vec,extra_var,case_study,...
        grad_step_size_vec,grad_dir_vec,num_grad_batches,grad_noise_lvl,i_b2b,...
        smooth_set_constr,U_out_vec,Obj_fun_vec);
    
    %Store ouptput measurements of past batches
    batch_output_mean_history(:,i_b2b) = proc_out_mean(:);
    batch_output_max_history(:,i_b2b) = proc_out_max(:);
    batch_output_min_history(:,i_b2b) = proc_out_min(:);
    
    %Store gradient measurement from past batches
    batch_gradient_history(i_b2b,:) = grad_noisy; 
    
%     disp(grad_additive);
%     disp(grad_noisy);
    %% %%%%%%%%%%%%%%%%% Parameter Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %number of samples for sampling uncertain parameter space
    num_samples = 1;
    
    %relative range for uncertain parameter space
    %(more stable to get the relative range instead of sampling joint
    %distributon)
    rel_range = 0.01;
    
    %number of selected parameter
    %(no treshold is used in this version)
    num_selected_paras = 2;
    
    %criteria for which sensitivities are considered in parameter ranking
    
    %consideration of output and gradient sensitivities
    sensitivity_crit = ['output','gradient'];

    %in case only output sensitivities are considered
%     sensitivity_crit = ['output'];
   
    
    %perform parameter selection
    para_selection = parameter_selection(U_vec,K_model,num_samples,...
    rel_range,num_selected_paras,case_study,sampling,num_outputs,sensitivity_crit,...
    Y_initial,i_b2b,extra_var,grad_dir_vec,grad_step_size_vec);
    
    %in case of fixed parameter subset (K_I and K_X) from Mandur and Budman (2015)
%     para_selection = [0 1 0 0 1 0 0 0];

    disp('Selected Parameters');
    disp(para_selection)


    %% %%%%%%%%%%%%%%%%% Step 2: Model Parameter Estimation %%%%%%%%%%%%%%%
    
    
    %introduce further limitation in the change in parameters. This is
    %necessary if a larger number of parameters is used for estimation as
    %parameters with less sensitivity lead to a more variability
    
    %parameter limitation in percentage from previous values 
    %(if para_limit = 0 then there is not extra limitation)
    %(if para_limit = 1 then parameters can change by 100%)
    %(if para_limit = .2 then parameters can change by 20%)
    para_limit = 1;
    
    %Indicate the number of past batches whose exprimental data are used for estimation
    num_batch_estim = 1;
    
    
    %call up function to fit model to measured process outputs
    [estim_model_para] = model_parameter_estimation(batch_output_mean_history,...
        proc_out_max,proc_out_min,U_vec,sampling,extra_var,K_model,...
        lb_model_parameter,ub_model_parameter,case_study,prev_corr,...
        para_selection,Y_initial,num_outputs,para_limit,num_batch_estim,i_b2b);

    %% %%%%%%%%%%%%%%%%% Step 3: Gradient Correction %%%%%%%%%%%%%%%%%%%%%%

    
    %indicate the number of past batches whose measurements are used in gradient correction    
    num_corrected_grad = 1;
    
    %for additional limitation in the parameters
    para_limit_grad = 0;
%     para_limit_grad = 0.75;

    %perform gradient correction
    [K_model_prime,new_corr,pred_grad] = gradient_correction(estim_model_para,batch_gradient_history,U_vec,i_b2b,...
    case_study,para_limit_grad,lb_model_parameter,ub_model_parameter,Y_initial,para_selection,...
    prev_corr,grad_dir_vec,grad_step_size_vec,sampling,extra_var,num_corrected_grad,eps_bound,...
    batch_output_max_history,batch_output_min_history,num_outputs);
     

% % in case of no gradient correction
%     new_corr = prev_corr;
%     K_model_prime = estim_model_para;
%     plant_grad_vec = 0;
%     grad_pred_vec = 0;
%     pred_grad = 0;


    %% %%%%%%%%%%%%%%%%% Step 4: Model-based Optimization %%%%%%%%%%%%%%%%%
    
    %perform a model-based optimization
    [U_new] = model_based_optimization(U_vec,K_model_prime,case_study,...
    lb_opt,ub_opt,reactor_vol_constr,extra_var,Y_initial,new_corr,sampling,i_b2b); 
    

    %update decision variable history vector
    U_vec(i_b2b + 1,:) = U_new;
    K_model = K_model_prime;
    K_model_vec(i_b2b,:) = K_model_prime; 
    
    %update correction
    prev_corr = new_corr;
    
    new_corr_mat(i_b2b,:) = new_corr(:);
  
    
      %% Plot Process and Model Outputs
    
            %plot in single plot
        plot_criterion = 'single';  
        %plot in separate plots
    %     plot_criterion = 'separate';

        %plot set-based constraints
        plot_set_based_constaints = 'on';

        plot_model_process_crit = {'process','model'};
        %Plot process output measurements
        plot_pen_process_outputs(proc_out_mean,sampling,proc_out_max,...
            proc_out_min,list_output,num_outputs,plot_criterion,...
            plot_set_based_constaints,plot_model_process_crit,estim_model_para,...
            Y_initial,U_vec,i_b2b,extra_var,prev_corr,new_corr)
        
        %plot predicted cost function
%         plot_pen_model_cost_function(pen_process_optimum,Y_initial,sampling,...
%         K_model_prime,extra_var,new_corr,U_vec,U_out_vec,Obj_fun_vec,i_b2b)
    

        %plot trajectory of decision variables
        plot_penicillin_decision_trajectory(U_vec,i_b2b,pen_process_optimum,[]);

    %% Output File
    %define output vector
    output_vec = [K_model_prime batch_gradient_history(i_b2b,:) pred_grad U_vec(i_b2b,:)];  

    %output file
    file_name = ['b2b_pen_para_selection_',num2str(i_sim),'.txt'];
    if i_b2b == 1
        dlmwrite(file_name, output_vec, 'delimiter', '\t','precision',  '%.4f');
        open(file_name);
    else
        dlmwrite(file_name, output_vec,'-append', 'delimiter', '\t','precision',  '%.4f');
    end
    
end
end

%% B2B Algorithm End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end