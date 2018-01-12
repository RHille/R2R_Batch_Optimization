function dY = penicillin_process_simulator(~, Y,U,K_sim,extra_var)
%mathematical equations to simulate the penicillin process

dY=zeros(4,1);

X=Y(1);
P=Y(2);
S=Y(3);
V=Y(4);

% Input Feed rate
F=U(2);

% Substrate conc. in the feed
sf=extra_var(1); % g/l

% Model Parameters
% K_sim=[0.092 0.15 0.005 0.0002 0.1 0.04 0.45 0.9 0.014]';

%To introduce process disturbances, noise can be introduced into the
%parameter values
% parameter_noise_lvl = extra_var(2)
% K_sim = K_sim.*(1+parameter_noise_lvl*randn(size(K_sim)));

% Biomass
mux=K_sim(1);
Kx=K_sim(2);
% Product 
mup=K_sim(3);
Kp=K_sim(4);
KI=K_sim(5);
Kh=K_sim(6);
% Substrate
Yxs=K_sim(7);
Yps=K_sim(8);
mx=K_sim(9);

% Loss in the culture volume due to evaporation
Floss=V*2.5*10^-4*(exp(5*(298-273)/(373-273))-1);

% ODEs
dY(1)=mux*S*X/(Kx*X+S)-(X/V)*(F-Floss);
dY(2)=mup*S*X/(Kp+S*(1+S/KI))-Kh*P-(P/V)*(F-Floss);
dY(3)=-mux*S*X/(Kx*X+S)/Yxs-mup*S*X/(Kp+S*(1+S/KI))/Yps-mx*X+F*sf/V-(S/V)*(F-Floss);
dY(4)=F-Floss;

end