function [model_check] = model_update_criterion(U,sampling,...
    Y_initial,para_model,y_out,prev_corr,y_max,y_min)


%define model used for process simulation
simple_model = @simple_process_model;


opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% process optimization

Y0 = Y_initial;

%initial parameter values
K_ini = para_model;

%get model prediction
[T, y_model] = ode45(simple_model,sampling,Y0,opt,para_model(1:2),U);
y_model = y_model - prev_corr;



myc = [y_model-y_max';y_min'-y_model];


if max(myc) >= 1e-2
    model_check = 0;
else
    model_check = 1;
end

end