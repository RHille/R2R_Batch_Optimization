function [f1,c1,ceq1,S,variance] = compute_robust_loop_input_uncertainty(U,U_nom,new_corr,K,theta1,theta2,epsi1,epsi2,w,pci_func,pci_den,para_selection,sampling,bio_nom,sf_nom_mod)

U_dist_diff = U(1) - U_nom(1);


nsamples=length(theta1);

%this computation requires a processor with four cores
spmd
    if (labindex==4)
        ai=(labindex-1)*floor(nsamples/numlabs)+1;
        if (rem(nsamples,numlabs)>0)
            bi=nsamples;
        else
            bi=labindex*floor(nsamples/numlabs);
        end
    else
        ai=(labindex-1)*floor(nsamples/numlabs)+1;
        bi=labindex*floor(nsamples/numlabs);
    end
end

% opt=odeset('NonNegative',[1,2,3,4],'RelTol',1e-6,'AbsTol',1e-8);
opt=odeset('RelTol',1e-6,'AbsTol',1e-8);

numofgridpts=length(theta1);
S=zeros(length(theta1),1);
constr=zeros(length(theta1),1);
Khat=K;
spmd
    for i=ai:bi
        
        %change selected parameter
%         i_loop_count = 1;
%         for i_loop1 = 1:length(K)
%             
%             if para_selection(i_loop1) == 1 && i_loop_count == 2 
%                 Khat(i_loop1) = theta2(i);
%             end
%             if para_selection(i_loop1) == 1 && i_loop_count == 1 
%                 Khat(i_loop1) = theta1(i);
%                 i_loop_count = i_loop_count + 1;
%             end
% 
%         end
%         Khat(2)=theta1(i);
%         Khat(5)=theta2(i);

                   U_dist(1) = double(theta1(i))+U_dist_diff;
                   bio_sample = double(theta2(i));
        %            sf_sample = double(theta2(i));
                   


        %         Khat(2)=theta1(i);
        %         Khat(5)=theta2(i);
                [~,dX] = ode15s(@process_opt_model_inhibit_input,sampling,[bio_sample 0 U_dist(1) 100],opt,[U_dist(1) U(2)],Khat,sf_nom_mod);

%         [~,dX] = ode15s(@process_opt_model_inhibit,sampling,[0.1 0 U(1) 100],opt,U,Khat);
        dX(:,:)=dX(:,:)-new_corr;
        S(i,:)=dX(end,2)*dX(end,4);
%         disp('S compute loop');
%         disp(S(i,:));
%         figure(93)
%         drawnow
%         scatter(S(i,:),i)
%         hold on;
%         grid on;
        constr(i,:)=dX(end,4);
    end
end

Shat=zeros(length(theta1),1);
chat=zeros(length(theta1),1);


for i=1:4
    Shat=Shat+S{i};
    chat=chat+constr{i};
end
S=Shat;
constr=chat;

% disp('S compute loop 2');
%         disp(S);
        
% figure(93)
% grid on;
% drawnow
% histogram(S,20,'Normalization','probability');
% for i_2 = 1:length(theta1)
%  scatter(S(i_2),i_2,'b');
% end
% hold off

a=(w'*0.25*S)/sum(w*0.25);
numofterms=length(pci_den)+1;

coeffs=(w'*(0.25*pci_func(epsi1,epsi2).*repmat(S,1,numofterms-1)))./pci_den;
error=sum((repmat(S,1,numofterms)-cumsum([repmat(a(1,1),numofgridpts,1) pci_func(epsi1,epsi2).*repmat(coeffs,numofgridpts,1)],2)).^2);
lim=find(error==min(error));
a=[a(1,1) coeffs(1:lim-1)];

% mean=a(1);
variance=sum((a(2:end).^2).*pci_den(1:lim-1));
% f1=-mean+0*variance;
%%%%%%%%%%%%
c1=max(constr)-120;
%%%%%%%%%%%
ceq1=[];

disp(U)

% opt=odeset('NonNegative',[1,2,3,4],'RelTol',1e-6,'AbsTol',1e-8);
opt=odeset('RelTol',1e-6,'AbsTol',1e-8);
% K(2)=Kout_prime(1); K(5)=Kout_prime(2); 
[~,dX] = ode15s(@process_opt_model_inhibit_input,sampling,[bio_nom 0 U(1) 100],opt,U,K,sf_nom_mod);
dX(:,:)=dX(:,:)-new_corr;
S=dX(end,2)*dX(end,4);
% constr=dX(end,4);
% ceq1=[];
% c1=constr-120;
% f1=-S+0.01*variance;
f1=-S+0.2*variance;
% f1=-S+0.05*variance;
% f1=-S;
end