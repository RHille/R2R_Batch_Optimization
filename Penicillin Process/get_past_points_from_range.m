function [grad_point_vec_mat]= get_past_points_from_range(K_model,initial_cond,U_vec,sampling,extra_var,...
    U_out_vec,i_b2b,num_grad_selec,num_grad_points,lb_doe,ub_doe,pen_std)



%define model
model = @penicillin_process_model;

Y0 = initial_cond;
Y0_grad1 = Y0;
Y0_grad2 = Y0;
U_curr = U_vec(i_b2b,:);
U_new = U_curr;
opt=odeset('RelTol',1e-6,'AbsTol',1e-8);



%run model for range of values
U_range = linspace(1,75,75);

%calculate nominal gradients
Y0_grad1(3) = U_curr(1);
[~, Y_nom] = ode15s(model,sampling, Y0_grad1, opt,U_curr,K_model,extra_var);

%objective funcation value at current operating point
obj_curr = Y_nom(end,2)*Y_nom(end,4);

for i = 1:length(U_range)
    
    U_new(1) = U_range(i);
    
    %calculate objective function value
    Y0_grad2(3) = U_new(1);
    [~, Y_grad] = ode15s(model,sampling, Y0_grad2, opt,U_new,K_model,extra_var);
    
    obj_val_vec(i) = Y_grad(end,2)*Y_grad(end,4);
    
end

%relative objective function vector
rel_obj_vec = obj_val_vec/obj_curr;

if i_b2b > 2
    %get range from measurement standard deviation
    rel_range = pen_std/obj_curr;
else
    %permissible relative range
    rel_range = 0.95;
end

%obtain range
cost_bound_counter = 1;
for i_cost_bound = 1:length(rel_obj_vec)
    if cost_bound_counter == 1 && rel_obj_vec(i_cost_bound) >= rel_range
        lb_range = U_range(i_cost_bound);
        cost_bound_counter = 2;
    elseif cost_bound_counter == 2 && rel_obj_vec(i_cost_bound) <= rel_range
        ub_range = U_range(i_cost_bound);
        cost_bound_counter = 3;
    end
end

%if doe is performed
if i_b2b > 2
    lb_range = lb_doe;
    ub_range = ub_doe;
end

grad_point_vec = [];

% grad_point_vec = U_out_vec >= lb_range;
% grad_point_vec = U_out_vec <= lb_range || U_out_vec >= ub_range;
% grad_point_vec = grad_point_vec.*(U_out_vec <= ub_range);
% grad_point_vec = grad_point_vec.*(U_out_vec >= ub_range);

grad_point_vec = zeros(size(U_out_vec,1),1);
for i_grad_p = 1:length(grad_point_vec)
    if U_out_vec(i_grad_p) <= lb_range+1 || U_out_vec(i_grad_p) >= ub_range-2
        grad_point_vec(i_grad_p) = 1;
    end
end

lin_range = fliplr(linspace(1,length(U_out_vec(:,1)),length(U_out_vec(:,1))));

grad_point_vec = grad_point_vec(:,1);

grad_point_vec = grad_point_vec.*(lin_range' <= num_grad_selec);

grad_point_vec_mat = [grad_point_vec U_out_vec];
end