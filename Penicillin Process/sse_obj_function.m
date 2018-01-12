function [K_model_new,SSE] = sse_obj_function(proc_out_mean_history,proc_out_max,proc_out_min,...
        U_vec,sampling,extra_var,K_model,lb,ub,case_study,prev_corr,para_selection,initial_cond,...
        num_outputs,para_limit,num_batch_estim,i_b2b)
%function to fit model parameters to measured process outputs

disp('PARAMETER ESTIMATION STEP');

%initial parameter values
K_ini = K_model;

%in case of further limitation in parameter values
if para_limit == 0
    lb_new = lb;
    ub_new = ub;
else
    rel_var = para_limit;
    lb_new=max(K_ini.*(1-rel_var),lb);
    ub_new=min(K_ini.*(1+rel_var),ub);
end

%initial process conditions
Y0 = initial_cond;

if strcmp(case_study,'penicillin')
    %Define penicillin process model
    model = @penicillin_process_model;
end


%ode options
opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% 'NonNegative',[1,2,3,4],


%fmincon options
options=optimset('Algorithm','interior-point','Display','iter','UseParallel','always',...
    'TolX',1e-4,'TolFun',1e-4,'TolCon',1e-4,'MaxFunEvals',1000);

% %'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-8
% [x,SSE] = fmincon(@sseobj,K_ini,[],[],[],[],lb_new,ub_new,[],options);

SSE = sseobj(K_ini);

K_model_new = K_ini;
% for i2 = 1:length(K_model)
%   if para_selection(i2) == 1
%      K_model_new(i2) = x(i2);
%   end
% end


function F = sseobj(x)
%objective function which calls up sse caluclation for each operating point
    
    %define overall sse
    sse_all = 0;
    
    %use data from the last num_batch_estim batches. however only if the
    %number of iterations i_b2b is equal or greater to num_batch_estim
    for i_sse_obj = 1:min(i_b2b,num_batch_estim)
        
      K_new = K_ini;
      %only change the parameters which have been selected for estimation
      for i_selection = 1:length(K_model)
          if para_selection(i_selection) == 1
              K_new(i_selection) = x(i_selection);
          end
      end
        
       
        %define operating point, start with the current one
        U_obj = U_vec(i_b2b+1-i_sse_obj,:);
        
        %get experimental data of corresponding batch
        proc_out_mean = reshape(proc_out_mean_history(:,i_b2b+1-i_sse_obj),length(sampling),num_outputs);
        
        %get sse for operating point U_obj
        sse_all(i_sse_obj) = sseobj_one_batch(K_new,U_obj);
        
    end

       F = sum(sse_all);
       
  function F2 = sseobj_one_batch(K_new,U_obj) 
      
      
        if strcmp(case_study,'penicillin')
        
        Y0(3) = U_obj(1);
        %update initial substrate concentration
        [T, Y] = ode45(model,sampling, Y0, opt,U_obj,K_new,extra_var);
    
        %add correction for accurate prediction
        Y=Y-prev_corr;

        %calculate the sum of square errors
        sse_vec = 0;
        
        W_sse = [4 0.4 4];
        
        for i_sse = 1:num_outputs-1
            
            %model fitting cost function
            sse = 1/abs(mean(proc_out_mean(:,i_sse)))*((Y(:,i_sse)-proc_out_mean(:,i_sse))).^2;
%             sse = 1/(W_sse(i_sse))*((Y(:,i_sse)-proc_out_mean(:,i_sse))).^2;
%             sse = ((Y(:,i_sse)-proc_out_mean(:,i_sse))).^2;
            sse_vec(i_sse) = sum(sse);
        end
        end  
       
        F2=sum(sse_vec);
  end
end


end