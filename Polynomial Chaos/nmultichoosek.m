function nmck = nmultichoosek(bins,balls) 
    n = bins;
    k = balls;
    nmck = nchoosek(n+k-1,n-1); 
end 