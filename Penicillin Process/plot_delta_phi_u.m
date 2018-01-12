function [P_final_process] = plot_delta_phi_u(U_vec,U_out_vec,Obj_fun_vec,i_b2b,num_grad_plot,initial_cond,...
    extra_var,model_para,sampling,plotting_mode,K_sim,P_final_process,grad_point_vec_out,grad_point_vec_mat)

min_plot = min(i_b2b,num_grad_plot)-1;

Y0 = initial_cond;
U = U_vec(i_b2b,:);

%define process simulator
model = @penicillin_process_model;
process = @penicillin_process_simulator;

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

if strcmp(plotting_mode,'points')

    S0_vec = U_out_vec(end-min_plot:end,1);
    for i = 1:length(S0_vec)

        U(1) = S0_vec(i);
        Y0(3) = U(1);
        %run model
        [~, Y] = ode45(model,sampling, Y0, opt,U,model_para,extra_var);
        P_final(i) = Y(end,2)*Y(end,4);
    end

    %plot
    figure(31)
    clf
    scatter(U_out_vec(end-min_plot:end,1)-U_out_vec(end,1),Obj_fun_vec(end-min_plot:end)-Obj_fun_vec(end),'fill');
    hold on
    scatter(U_out_vec(end-min_plot:end,1)-U_out_vec(end,1),P_final(end-min_plot:end)-P_final(end),'r','fill');
    grid on;
    drawnow
    hold off

elseif strcmp(plotting_mode,'all') && i_b2b > 1
    
    S0_vec = linspace(0.1,70,30);
    for i = 1:length(S0_vec)

        U(1) = S0_vec(i);
        Y0(3) = U(1);
        %run model
        [~, Y1] = ode45(model,sampling, Y0, opt,U,model_para,extra_var);
        P_final_model(i) = Y1(end,2)*Y1(end,4);
        %run process
        
%         if i_b2b == 1
            [~, Y2] = ode45(process,sampling, Y0, opt,U,K_sim,extra_var);
            P_final_process(i) = Y2(end,2)*Y2(end,4);
%         end
    end
    
    U(1) = U_out_vec(end,1);
    Y0(3) = U(1);
    [~, Y1] = ode45(model,sampling, Y0, opt,U,model_para,extra_var);
    P_final_model_curr = Y1(end,2)*Y1(end,4);
    [~, Y2] = ode45(process,sampling, Y0, opt,U,K_sim,extra_var);
    P_final_process_curr = Y2(end,2)*Y2(end,4);
    
    %plot
    figure(31)
    clf
    plot(S0_vec,P_final_model-P_final_model_curr,'b')
    
    hold on
    plot(S0_vec,P_final_process-P_final_process_curr,'r')
    scatter(U_out_vec(end-min_plot:end,1),Obj_fun_vec(end-min_plot:end)-Obj_fun_vec(end),'r','fill');
    
    %remove gradient measurements
    
    grad_point_vec_out = [grad_point_vec_out(1:end-2,:); 1];
    
%     size(grad_point_vec_out)
%     size(U_out_vec)
% %     for i_u_out = 1:length(grad_point_vec_out)
% %         if grad_point_vec_out(end+1-i_u_out) == 1
% %             scatter(U_out_vec(end+1-i_u_out,1),Obj_fun_vec(end+1-i_u_out)-Obj_fun_vec(end),'b','fill');
% %         end
% %     end

 for i_u_out = 1:length(grad_point_vec_mat(:,1))
        if grad_point_vec_mat(end+1-i_u_out,1) == 1
            scatter(grad_point_vec_mat(end+1-i_u_out,2),Obj_fun_vec(end+1-i_u_out)-Obj_fun_vec(end),'b','fill');
        end
 end
  
%     scatter(U_out_vec(end-min_plot:end,1),Obj_fun_vec(end-min_plot:end)-Obj_fun_vec(end),'r','fill');
%     scatter(U_out_vec(end-min_plot:end,1)-U_out_vec(end,1),P_final(end-min_plot:end)-P_final(end),'r','fill');
    grid on;
    drawnow
    hold off
    
end


end