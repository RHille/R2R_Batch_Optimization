function model_adequacy(U,sampling,Y_initial,c_model,process_optimum)

%function to calculate the parameter values for which the model is
%adequate (predict an optimum at the process optimum)

%define model used for process simulation
simple_model = @simple_process_model;

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);



Y0 = Y_initial;

%generate vectors of parameter values
a_vec = [0.1:0.1:5];
b_vec = [0.1:0.1:5];


%input perturbation for gradient calculation
perturbation = 0.1;

%optimal input
U = process_optimum;

%run through loops to calculate the gradient
for i = 1:length(a_vec)
    for i2 = 1:length(b_vec)
        
        %calculate gradient at optimal input
        K_model = [a_vec(i) b_vec(i2)];
        grad_opt_input(i,i2) = get_gradient(K_model);
        
    end
end

%only consider gradients close to 0
% grad_opt_input = grad_opt_input.*(abs(grad_opt_input)<=0.1);
grad_opt_input = (abs(grad_opt_input)<=0.05);

[X,Y] = meshgrid(b_vec,a_vec);
%create contout plot
figure(4)
contour(b_vec,a_vec,abs(grad_opt_input),1)
title('model adequacy region');
% surf(X,Y,abs(grad_opt_input))


    function F = get_gradient(x)
        
        %calculate predicted gradient
        
        %nominal
        [T, y_model_nominal] = ode45(simple_model,sampling,Y0,opt,x,U);
        
        U_perturbed = U + perturbation;
        [T, y_model_perturbed] = ode45(simple_model,sampling,Y0,opt,x,U_perturbed);
        

        %calculate nominal cost
        cost_nominal = y_model_nominal(end) - c_model*U^2;

        %perturbed cost
        cost_perturbed = y_model_perturbed(end) - c_model*U_perturbed^2;

        %estimate gradient
        gradient_model = (cost_perturbed - cost_nominal)/perturbation;
        
        F = gradient_model;
    end

save('model_ade_cont.mat','a_vec','b_vec','grad_opt_input')
end