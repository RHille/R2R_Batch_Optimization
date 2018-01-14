function dY = simple_process_model(t, Y,Para_model,Y_in)
%mathematical equations to simulate the penicillin process

dY=zeros(1,1);

%define simple process
dY(1) = Para_model(1)*Y_in - Para_model(2)*Y(1); 

% dY(1) = Para_sim(1)*Y_in - Para_sim(2)*Y(1)^2 + Para_sim(3)*Y(1); 



end