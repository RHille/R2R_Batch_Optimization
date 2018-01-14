function B2B_syn_batch_process_doe
%Simple example to illustrate the algorithm for simultaneous identification
%and optimization using set-based constraints and an extended gradient
%correction with design of experiments

%Thesis, chapter 5, synthetic batch process

%Author: Rubin Hille


%Program Initialization 
clc
close all
clear

%% Process and Model Definition

%process parameters
para_sim = [1 2 0.3];

%initial model parameters
para_model_ini = [4 2 0.05];

%cost function parameter
cost_para = 0.05;
cost_para_model = 0.05;

%sampling time
sampling = [0:.5:5];
% sampling2 = [0:.1:5];

%initial value of decision variable
U0 = 20;

%initial batch condition
Y_initial = 8;

%define model
simple_model = @simple_process_model;

%% Simple Process Optimization and Illustration

%For test purposes, find optimum and plot objective function for the "true" process
%call up function find optimum
process_optimum = optimize_simple_process(U0,sampling,Y_initial,para_sim,cost_para)


%% Model Adequacy Check

% compute parameter space for which the model is adequate (predict the
% optimum at the correction location)
%  model_adequacy(process_optimum,sampling,Y_initial,cost_para_model,process_optimum)

%% B2B Algorithm Options

%amount of measurement noise (betweeen 0 and 1)
noise_lvl = 1;

%TYPE OF GRADIENT CORRECTION METHOD
%OPTION OF STANDARD OR EXTENDED GRADIENT CORRECTION

%in case of only local gradient correction with fixed input perturbations
% grad_corr_method = 'standard';

%in case of extended gradient correction with perturbations based DoE
grad_corr_method = 'extended';


%% B2B Algorithm Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i_sim = 1:10

%initialize values
Y0 = Y_initial;
U = U0;
U_vec = U;

%number of batch runs 
num_iterations = 50;

%initial correction term
prev_corr = 0;

%initial input perturbation for gradient measurements (will be updated if
%experimental design procedure is used)
perturbation = 1.5;

%bound on the relative truncation error
trunc_error = 0.01;

%parameters
para_vec = para_model_ini;
para_model = para_model_ini;

%initial cost and input matrices
U_out_vec = [];
Obj_fun_vec = [];

%initial bounds for experimental design (will be updated if
%experimental design procedure is used)
lb_doe = U;
ub_doe = U;


for i_b2b = 1:num_iterations

%% %%%%%%%%%%%%%%%%% Step 1: Process simulation %%%%%%%%%%%%%%%%%%%%%%%

%run process simulator
[output,gradient,y_max_sf,y_min_sf,U_out_vec,Obj_fun_vec,cost_std] = simple_process_simulation(para_sim,sampling,Y_initial,...
    U,noise_lvl,perturbation,cost_para,U_out_vec,Obj_fun_vec);


[U_out_vec Obj_fun_vec]
%% %%%%%%%%%%%%%%%%% Step 2: Parameter Estimation %%%%%%%%%%%%%%%%%%%%%%%

%model update criterion
[model_check] = model_update_criterion(U,sampling,Y_initial,para_model,output,...
    prev_corr,y_max_sf,y_min_sf)
% 
if model_check == 0
% if i_b2b == 1
%perform parameter estimation
[para_new] = model_parameter_estimation(U,sampling,Y_initial,para_model,output,...
    prev_corr,y_max_sf,y_min_sf);
end

%update parameters
para_model = para_new;

%plot process and predicted outputs
figure(1)
plot(sampling,output,'bo')
hold on
%simulate process
[T, y_model] = ode45(simple_model,sampling,Y0,[],para_new(1:2),U);
y_model = y_model - prev_corr;
plot(sampling,y_model,'r-','linewidth',1.8)
plot(sampling,y_max_sf,'--m','linewidth',1.8)
plot(sampling,y_min_sf,'--m','linewidth',1.8)
title('Model output fitting with sb const');
grid on;
xlabel('Time')
ylabel('Y')
hold off
drawnow


%% %%%%%%%%%%%%%%%%% Step 2: Gradient Correction %%%%%%%%%%%%%%%%%%%%%%%

%for different gradient correction methods
if ismember('standard',grad_corr_method)
    %standard (only local) gradient correction
    num_grad_corr = 2;
else
    %extended gradient correction
    num_grad_corr = min(5,length(U_out_vec));
end


grad_point_vec = get_grad_point_vec(U_out_vec,num_grad_corr,lb_doe,ub_doe)

%perform gradient correction
[para_corr,new_corr,pred_gradient,plant_gradient2] = gradient_correction(U,sampling,...
    Y_initial,para_model,gradient,prev_corr,perturbation,cost_para,cost_para_model,trunc_error,...
    U_out_vec,Obj_fun_vec,num_grad_corr,i_b2b,grad_point_vec);

%update parameter
para_new = para_corr;

%update parameters
para_model = para_new;

% % %w/o gradient correction
% % new_corr = prev_corr;
% % pred_gradient = 0;

%plot process and predicted cost function
plot_cost_function(U0,sampling,Y_initial,para_sim,cost_para,para_new,cost_para_model,new_corr,U_out_vec,Obj_fun_vec)
%% %%%%%%%%%%%%%%%%% Step 3: Process Optimization %%%%%%%%%%%%%%%%%%%%%%%

%perform model-based optimization
[U_new] = model_based_optimization(U,sampling,Y0,para_new,cost_para_model);

U_vec(i_b2b+1) = U_new;
U = U_new;
prev_corr = new_corr;

%plot decision variable trajectory
iter_vec = 0:i_b2b;
optim_vec = process_optimum*ones(i_b2b+1,1);
figure(3)
plot(iter_vec,optim_vec,'k--','linewidth',1.8)
hold on
plot(iter_vec,U_vec,'r-o','linewidth',1.8)
title('Decision variable trajectory');
legend('actual optimal input','predicted optimal input');
grid on;
% hold off
drawnow

para_vec(i_b2b+1,:) = para_new;


%% %%%%%%%%%%%%%%%%% Step 4: Design of Experiments %%%%%%%%%%%%%%%%%%%%%%%
%for different gradient correction methods
if ismember('standard',grad_corr_method)
    %standard (only local) gradient correction
    doe_out = perturbation;
else
    %extended gradient correction
    [doe_out,lb_doe,ub_doe] = experimental_design_simple_example(U,sampling,Y_initial,para_model,...
    grad_point_vec,U_out_vec,num_grad_corr,cost_std);

    %update input perturbation
    perturbation = doe_out;
end

%% Output File
    %define output vector 
    output_vec = [para_new U_new plant_gradient2 pred_gradient model_check doe_out cost_std];  
    %output file
    file_name = ['syn_batch_process_doe_',num2str(i_sim),'.txt'];
%     file_name = ['simple_example_output_wo_doe_.txt'];
    if i_b2b == 1
        
        dlmwrite(file_name, output_vec, 'delimiter', '\t','precision',  '%.4f');
        open(file_name);
    else
        dlmwrite(file_name, output_vec,'-append', 'delimiter', '\t','precision',  '%.4f');
    end

end
end
%% B2B Algorithm End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save('optim_results_sio.mat','U_vec','iter_vec','optim_vec','para_vec')


end