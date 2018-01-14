function [output,gradient,y_max_sf,y_min_sf,U_out_vec,Obj_fun_vec,cost_std] = simple_process_simulation(K_sim,sampling,...
Y_initial,U,noise_lvl,perturbation,c,U_out_vec,Obj_fun_vec)

%define model used for process simulation
simple_simulator = @simple_process_simulator;

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

%get initial value
Y0 = Y_initial;

%inlet concentration based on decision variable
Y_in = U;

%run numerical ode solver 
[T, Y] = ode45(simple_simulator,sampling,Y0,opt,K_sim,Y_in);

%standard deviation of measurement noise
std_noise = 0.5*noise_lvl;

num_run = 4;

for i_run = 1:num_run
    %add measurement noise
    Y_noisy(i_run,:) = normrnd(Y,std_noise,size(Y,1),1);
        
end

%define output
% output = Y_noisy(1,:)';
output = mean(Y_noisy(:,:))';

%get standard deviation in cost measurement
cost_noisy = Y_noisy(:,end) - c*Y_in^2;

cost_std = sqrt(var(cost_noisy));

%% just in case of set-based constraints


%add smoothing factor
[y_max_sf,y_min_sf] = smoothing_factor(sampling,Y_noisy);


%% estimate gradient

%input perturbation
Y_in_perturb = Y_in + perturbation;

%run numerical ode solver 
[T, Y_perturbed] = ode45(simple_simulator,sampling,Y0,opt,K_sim,Y_in_perturb);

%add measurement noise
Y_perturbed_noisy = normrnd(Y_perturbed,std_noise,size(Y,1),1);

%calculate nominal cost
cost_nominal = output(end) - c*Y_in^2;
% cost_nominal = Y_perturbed_noisy(end) - c*Y_in^2;

%perturbed cost
cost_perturbed = Y_perturbed_noisy(end) - c*Y_in_perturb^2;

%estimate gradient
gradient = (cost_perturbed - cost_nominal)/perturbation;

U_out_vec = [U_out_vec; Y_in_perturb; Y_in];
Obj_fun_vec = [Obj_fun_vec; cost_perturbed; cost_nominal];


    function [y_max_sf,y_min_sf] = smoothing_factor(sampling,y_noisy_mat)
        
        %calculate a simple smoothing factor to add to the mean
        y_max_factor(1) = 0;
        y_min_factor(1) = 0;
        y_mean(1) = mean(y_noisy_mat(:,1));
        y_max_nom_curr(1) = max(y_noisy_mat(:,1));
        y_min_nom_curr(1) = min(y_noisy_mat(:,1));
        
        for i=2:length(sampling)
            
            %calculate mean
            y_mean(i) = mean(y_noisy_mat(:,i));
            y_max_nom_curr(i) = max(y_noisy_mat(:,i));
            y_min_nom_curr(i) = min(y_noisy_mat(:,i));
            
            %difference of mean to max and min vals
            y_max_factor(i) = max(y_max_factor(i-1),y_max_nom_curr(i)-y_mean(i));
            y_min_factor(i) = max(y_min_factor(i-1),y_mean(i)-y_min_nom_curr(i));  
        end
        
        y_max_sf = y_mean+max(y_max_factor);
        y_min_sf = y_mean-max(y_min_factor);
        
    end
end