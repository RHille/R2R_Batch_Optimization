function dY = penicillin_sens_analysis_model(t,Y,flag,K,U,extra_var)
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


% % Biomass
mux=K(1);
Kx=K(2);
% Product 
mup=K(3);
Kp=K(4);
KI=K(5);

% Substrate
Yxs=K(6);
Yps=K(7);
mx=K(8);

% Loss in the culture volume due to evaporation
Floss=V*2.5*10^-4*(exp(5*(298-273)/(373-273))-1);

% ODEs
dY(1)=mux*S*X/(Kx*X+S)-(X/V)*(F-Floss);
dY(2)=mup*S*X/(Kp+S*(1+S/KI))-(P/V)*(F-Floss);
dY(3)=-mux*S*X/(Kx*X+S)/Yxs-mup*S*X/(Kp+S*(1+S/KI))/Yps-mx*X+F*sf/V-(S/V)*(F-Floss);
dY(4)=F-Floss;

end