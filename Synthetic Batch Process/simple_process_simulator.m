function dY = simple_process_simulator(t, Y,Para_sim,Y_in)
%mathematical equations to simulate the penicillin process

dY=zeros(1,1);

%define simple process
dY(1) = Para_sim(1)*Y_in - Para_sim(2)*Y(1) + Para_sim(3)*t; 

% dY(1) = Para_sim(1)*Y_in - Para_sim(2)*Y(1)^2 + Para_sim(3)*Y(1); 



end