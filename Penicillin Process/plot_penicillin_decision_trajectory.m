function plot_penicillin_decision_trajectory(U_vec,i_b2b,process_opt,U_doe_vector)
%function to plot the trajectory of decision variables for the penicillin
%process

%Due to the constraint on the volume of the reactor, convergence to the
%optimal flow rate is achieved within one iteration. However, it takes
%several batch runs for the initial substrate concentration to converge to
%the optimum. For that reason, only the trajectory of the initial substrate
%concentration will be plotted

figure(2)
hold on;
grid on;
plot(1:i_b2b+1,process_opt(1)*ones(i_b2b+1,1),'--k','LineWidth',1.25);
plot(1:i_b2b+1,U_vec(:,1),'-kd','LineWidth',1.25);

if i_b2b > 1
    scatter(2:length(U_doe_vector)+1,U_doe_vector,'r','*')
end
ylabel('S_0');
xlabel('Iteration');
drawnow


end