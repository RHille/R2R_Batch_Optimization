function ind = ndind(x) 
% NDIND Compute the index of an  
 
if (~isvector(x)) 
    ind = zeros(1,size(x,2)); 
    for (kk=1:size(x,2)) 
        ind(kk) = ndind(x(:,kk)); 
    end 
    return 
end 

if (length(x) == 1) 
    ind = x(1); 
else 
    level = norm(x,1); 
    count = 0;
    for i=0:level-1
        count = count + nmultichoosek(length(x), i);
    end
    ind = count + ndind(x(2:end)); 
end 
     
end