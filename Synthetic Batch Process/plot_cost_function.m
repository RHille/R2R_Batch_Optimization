function plot_cost_function(initial_decision_val,sampling,...
    Y_initial,Para_sim,cost_para,para_model,c_model,prev_corr,U_out_vec,Obj_fun_vec)


%intial values for initial concentration
U = initial_decision_val;

%define model used for process simulation
simple_simulator = @simple_process_simulator;
simple_model = @simple_process_model;


opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% process optimization

Y0 = Y_initial;

%parameter of cost function
c = cost_para;

%define vector for Y_in
Y_in_vec = [0.1:0.5:15];

for i=1:length(Y_in_vec)
    
    %evaluate objective function
    obj_val(i) = objfun(Y_in_vec(i));
    obj_val_model(i) = objfun_model(Y_in_vec(i));
    
end

%plot objective function
figure(2)
plot(Y_in_vec,obj_val,'b','LineWidth',2)
hold on
plot(Y_in_vec,obj_val_model,'r','LineWidth',2)
scatter(U_out_vec,-Obj_fun_vec,'b','filled');
title('Cost function prediction');
legend('actual cost','predicted','measured');
grid on;
xlabel('Y_{in}')
ylabel('Cost')
hold off

%save cost function data
% save('cost_function_data_sio.mat','Y_in_vec','obj_val','obj_val_model')

    function F = objfun(x)
        
        %run numerical ode solver 
        [T, Y] = ode45(simple_simulator,sampling,Y0,opt,Para_sim,x);
        
        %cost function
        F = -(Y(end) - c*x^2);
        
    end

   function F = objfun_model(x)
        
        %run numerical ode solver 
        [T, Y] = ode45(simple_model,sampling,Y0,opt,para_model(1:2),x);
        Y = Y - prev_corr;
        
        %cost function
        F = -(Y(end) - para_model(3)*x^2);
        
    end

end