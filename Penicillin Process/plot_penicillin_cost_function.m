function plot_penicillin_cost_function(process_optimum,initial_cond,...
    initial_decision_val,extra_var,sampling,K_process)
%function to create a plot of the actual penicillin process cost function

sf_nom = extra_var(1);

%define process simulator
pen_simulator = @penicillin_process_simulator;

%define grid sizes for plotting
S0_step = (process_optimum(1) - initial_decision_val(1))/10;
F_step = (process_optimum(2) - initial_decision_val(2))/10;

S0_vec = [initial_decision_val(1):S0_step:process_optimum(1)+3*S0_step];
F_vec = [initial_decision_val(2):F_step:process_optimum(2)+3*F_step];

disp('GENERATE PENICILLIN PROCESS COST FUNCTION');

opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

Y0 = initial_cond;
    i3 = 1;
    for i1 = S0_vec
        i4 = 1;
        for i2 = F_vec


            U = [i1 i2];
            Y0(3) = U(1);
            
            [~, Y] = ode15s(pen_simulator,[sampling], Y0, opt,U,K_process,sf_nom);
            
            
            S=Y(end,2)*Y(end,4);

            P_final_mat(i3,i4) = S;

            i4 = i4 + 1;
        end
        disp('.');
        i3 = i3 + 1;
    end
    
    %surface plot
    figure(1)
    surf(F_vec,S0_vec,P_final_mat)
    xlabel('F');
    ylabel('S0');
    zlabel('P(t_f)');
    title('Penicillin Cost Function Surface Plot')
    
    %contour plot
    figure(2)
    contour(F_vec,S0_vec,P_final_mat,20)
    xlabel('F');
    ylabel('S0');
    zlabel('P(t_f)');
    title('Penicillin Cost Function Contour Plot')


end