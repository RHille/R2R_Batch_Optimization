function plot_pen_model_cost_function(process_optimum,initial_cond,sampling,...
    model_para,extra_var,new_corr,U_vec,U_out_vec,Obj_fun_vec,i_b2b)
%function to create a plot of the penicillin process cost function


%define process simulator
model = @penicillin_process_model;

U = U_vec(i_b2b,:);

%define grid sizes for plotting
S0_step = 2;
S0_vec = [1:S0_step:70];

% tf_vec = [1:1:50];
% Gln_vec = tf_vec;
disp('GENERATE CHO PROCESS COST FUNCTION');

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

Y0 = initial_cond;
    i3 = 1;
    for i1 = S0_vec
%         i4 = 1;
%         for i2 = Gln_vec


%             U = [i1 i2];
            %update glucose concentration
            Y0(3) = i1;
            U(1) = i1;
%             Y0(5) = i2;
            
            [~, Y] = ode45(model,sampling, Y0, opt,U,model_para,extra_var);
            
%             F2=Y(end,2)*Y(end,4);
            Y = Y - new_corr;
            
            F=Y(end,2)*Y(end,4);
            
%             F= interp1(sampling,Y(:,11),i2)/i2;
            
%             F = sum(Y(:,11)./Y(:,2));
            
%             F = Y(end,11)/mean(Y(:,2)+1);
            

%             Mab_final_mat(i3,i4) = F;
            Pen_final_mat(i3) = F;
%             Pen_final_mat2(i3) = F2;

%             i4 = i4 + 1;
%         end
        disp('.');
        i3 = i3 + 1;
    end
    
      %surface plot
    figure(10)
    plot(S0_vec,Pen_final_mat)
    hold on
    scatter(U_out_vec(:,1),Obj_fun_vec,'fill')
    hold off
    xlabel('S0');
    ylabel('P(tf)');
    title('Pen Cost Function')
    
%     %surface plot
%     figure(13)
%     surf(Gln_vec,Glc_vec,Mab_final_mat)
%     xlabel('Gln');
%     ylabel('Glc');
%     zlabel('MAb');
%     title('MAb Cost Function Plot')
%     
%     %contour plot
%     figure(12)
%     contour(Gln_vec,Glc_vec,Mab_final_mat,20)
%     xlabel('Gln');
%     ylabel('Glc');
%     zlabel('MAb');
    title('MAb Cost Function Contour Plot')
end