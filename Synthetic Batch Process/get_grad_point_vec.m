function grad_point_vec = get_grad_point_vec(U_out_vec,num_grad_corr,lb_doe,ub_doe)

grad_point_vec = zeros(length(U_out_vec),1);

i_count = 1;
for i=1:length(U_out_vec)
    
    if i_count <= num_grad_corr
        if U_out_vec(end-(i-1)) <= lb_doe || U_out_vec(end-(i-1)) >= ub_doe
            grad_point_vec(end-(i-1)) = 1;
            i_count = i_count +1;
        end
    else
        break;
    end

    
end

end