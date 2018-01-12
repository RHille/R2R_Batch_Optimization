function plot_pen_process_outputs(process_output,sampling,process_max,...
    process_min,output_list,num_output,plot_crit,plot_set_based_constr,...
    plot_model_process_crit,K_model,Y_initial,U_vec,i_b2b,extra_var,...
    prev_corr,new_corr)
%function to process measurements

%simulate model (in case it is required)
if ismember('model',plot_model_process_crit);
    %simulate model
    model = @penicillin_process_model;
    Y0 = Y_initial;
    U = U_vec(i_b2b,:);
    Y0(3) = U(1);
    opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
    %update initial substrate concentration
    [T, model_output] = ode45(model,sampling, Y0, opt,U,K_model,extra_var);
    
    %add correction for accurate prediction
    model_output_prev=model_output-prev_corr;
    
    model_output_new=model_output-new_corr;
end

%plot handle for plotting legend
plot_handle_all = [];

if strcmp(plot_crit,'separate');
    for i_plot = 1:3

        figure(i_plot)
    %     hold on
        plot(sampling,process_output(:,i_plot),'LineWidth',1.5);
        grid on;
        legend(output_list(i_plot));
    end
elseif strcmp(plot_crit,'single');
        figure(1)
        clf
        hold on
        
    for i_plot = 1:3

       
        if ismember('process',plot_model_process_crit);
            plot_handle = plot(sampling,process_output(:,i_plot),'LineWidth',1.5);
            if strcmp(plot_set_based_constr,'on')
%                 plot(sampling,process_max(:,i_plot).*(process_max(:,i_plot)>0),'--r','LineWidth',1);
%                 plot(sampling,process_min(:,i_plot).*(process_min(:,i_plot)>0),'--r','LineWidth',1);
                plot(sampling,process_max(:,i_plot),'--r','LineWidth',1);
                plot(sampling,process_min(:,i_plot),'--r','LineWidth',1);
            end
        end
        if ismember('model',plot_model_process_crit);
%             plot(sampling,model_output_prev(:,i_plot),'-.k','LineWidth',1);
            plot(sampling,model_output_new(:,i_plot),'-.b','LineWidth',1);
%             plot(sampling,model_output(:,i_plot),'-.k','LineWidth',1);
        end
        
        %update plot handle
        plot_handle_all = [plot_handle_all,plot_handle];
    end
        
        hold off
        grid on;
        legend([plot_handle_all],output_list);
        drawnow

end


save('pen_sb_info.mat','process_max','process_min','sampling','model_output_new');




end

