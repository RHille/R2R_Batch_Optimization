function [theta1_PCE, theta2_PCE,epsi_1,epsi_2,F] = twoP_legendre_sse_normal_shiftedmeans(Kout_prime,cov_data,pci_epsi1,pci_epsi2,num_poly)
% clear all; clc; 
%close all
% syms epsi
% for i=1:2
%     epsi(i)=sym(strcat('epsi',num2str(i)),'real');
% end

% epsis = [sym('epsi1','real') sym('epsi2','real')];
% 
% pci_epsi1=pcepoly1d('legendre',epsis(1),20);
% pci_epsi2=pcepoly1d('legendre',epsis(2),20);

% j=1;
% z=1;
% output=dlmread(['outputdata_6mar_RBTB_trunc5_noise_exp',num2str(z),'.txt']);
% 
% mu = output(j,5:6);
% Sigma = ([output(j,7) output(j,9); output(j,9) output(j,8)]);

K_data=[Kout_prime cov_data];

gleq=[1 	-0.997263861849 	0.00701814576495
2 	-0.985611511545 	0.0162774265831
3 	-0.964762255588 	0.0253910098329
4 	-0.934906075938 	0.0342745478477
5 	-0.896321155766 	0.0428359896785
6 	-0.849367613733 	0.0509978738117
7 	-0.794483795968 	0.0586839394615
8 	-0.73218211874 	0.0658220603578
9 	-0.66304426693 	0.0723456094297
10 	-0.587715757241 	0.078193695762
11 	-0.506899908932 	0.083311711103
12 	-0.421351276131 	0.0876518688047
13 	-0.331868602282 	0.0911736454878
14 	-0.239287362252 	0.0938441590423
15 	-0.144471961583 	0.0956384754512
16 	-0.0483076656877 	0.0965398415811
17 	0.0483076656877 	0.0965398415811
18 	0.144471961583 	0.0956384754512
19 	0.239287362252 	0.0938441590423
20 	0.331868602282 	0.0911736454878
21 	0.421351276131 	0.0876518688047
22 	0.506899908932 	0.083311711103
23 	0.587715757241 	0.078193695762
24 	0.66304426693 	0.0723456094297
25 	0.73218211874 	0.0658220603578
26 	0.794483795968 	0.0586839394615
27 	0.849367613733 	0.0509978738117
28 	0.896321155766 	0.0428359896785
29 	0.934906075938 	0.0342745478477
30 	0.964762255588 	0.0253910098329
31 	0.985611511545 	0.0162774265831
32 	0.997263861849 	0.00701814576495];

mu = [K_data(1) K_data(2)];
Sigma = [K_data(3) K_data(5); K_data(5) K_data(4)];

chi2val=chi2inv(0.999,2);
rect_coord=[mu(1)-sqrt(chi2val*Sigma(1,1)), mu(2)-sqrt(chi2val*Sigma(2,2)), 2*sqrt(chi2val*Sigma(1,1)), 2*sqrt(chi2val*Sigma(2,2))];
x_1=rect_coord(1)+rect_coord(3);
x_2=rect_coord(2)+rect_coord(4);

if rect_coord(1)<0
    rect_coord(1)=0;
end
if rect_coord(2)<0
    rect_coord(2)=0;
end
epsi_1 = linspace(rect_coord(1), x_1,100);
epsi_2 = linspace(rect_coord(2), x_2,100);

[X1,X2] = meshgrid(epsi_1,epsi_2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(epsi_2),length(epsi_1));



% figure(90)
% surf(X1,X2,F);
% drawnow


% figure(91)
% contour(epsi_1,epsi_2,F);
% % contour(F)
% grid on;
% drawnow
% hold on


theta1=X1; theta2=X2; p_theta=F;

mp_theta1=[];
for i=1:length(epsi_1)
    mp_theta1(i)=trapz(theta2(:,i),p_theta(:,i));
end

cum_p=trapz(theta1(1,:),mp_theta1);
p_theta=p_theta/cum_p;

mp_theta1=[];
for i=1:length(epsi_1)
    mp_theta1(i)=trapz(theta2(:,i),p_theta(:,i));
end

cp_theta2=p_theta./repmat(mp_theta1,size(p_theta,1),1);

cump_theta1=[]; cump_theta1(1)=0;
for j=2:length(mp_theta1)
    cump_theta1(j)=trapz(theta1(1,1:j),mp_theta1(1:j));
end

epsi1=gleq(:,2);
% theta1=interp1(cump_theta1',theta1(1,:)',unifcdf(epsi1,-1,1),'linear','extrap');

theta1_i=theta1(1,cump_theta1<=0.9999);
cump_theta1_i=cump_theta1(cump_theta1<=0.9999);
theta1=interp1(cump_theta1_i',theta1_i',unifcdf(epsi1,-1,1),'linear','extrap');


% cump_new = interp1(theta1,unifcdf(epsi1,-1,1),theta1_i','v5cubic');
cump_new = interp1(theta1,unifcdf(epsi1,-1,1),theta1_i','linear','extrap');

% pci=pcepoly1d('legendre',epsi(1),20);
pci=pci_epsi1;

% load c.mat
% chat=c;
c=[]; err=[];
c(1,1)=(gleq(:,3)'*0.5*theta1)/sum(gleq(:,3)*0.5);
err(1,1)=sum((theta1-repmat(c(1,1),length(theta1),1)).^2);

for i=2:num_poly
c(i,1)=sum(gleq(:,3).*subs(pci(i)).*theta1*0.5)/sum(gleq(:,3).*subs(pci(i)^2)*0.5);
err(i,1)=sum((theta1-repmat(c(1,1),length(theta1),1)-subs(pci(2:i))*c(2:i,1)).^2);

% if(err(i,:)<1e-6)
%     pci=pci(1:i); 
%     break;
% elseif(err(i,:)-err(i-1,:)>1)
%     pci=pci(1:i-1);
%     c=c(1:i-1,1);
%     err=err(1:i-1,1);
%     break;
% end
end
% disp(i)

% lim=find(err==min(err));
% c=c(1:lim,:);
% pci=pci(1:lim);
% max(max(c-chat))
% min(min(c))
theta1=pci*c; % theta1=c;
% save c.mat c


% cump_new_pci = interp1(double(subs(theta1)),unifcdf(epsi1,-1,1),theta1_i','v5cubic','extrap');
cump_new_pci = interp1(double(subs(theta1)),unifcdf(epsi1,-1,1),theta1_i','linear','extrap');
% cump_new_old = interp1(theta1_old,unifcdf(epsi1,-1,1),theta1_i','pchip','extrap');

figure(30)
drawnow
title('PCE Approximation');
grid on
plot(theta1_i(1:end-1),diff(cump_new_pci));
hold on
plot(theta1_i(1:end-1),diff(cump_new));
% plot(theta1_i(1:end-1),diff(cump_theta1_i));
hold off




cump_theta2=zeros(size(p_theta,1),size(p_theta,2));

for i=1:size(cp_theta2,2)
    cump_theta2(1,i)=0;
for j=2:size(cp_theta2,1)
    cump_theta2(j,i)=trapz(theta2(1:j,i),cp_theta2(1:j,i));
end
end

epsi2=gleq(:,2);
theta2_new=[];

for i=1:size(cp_theta2,2)
    theta2_i=theta2(cump_theta2(:,i)<=0.999,i);
    cump_theta2_i=cump_theta2(cump_theta2(:,i)<=0.999,i);
    theta2_new(:,i)=interp1(cump_theta2_i,theta2_i,unifcdf(epsi2,-1,1),'linear','extrap');
end
theta2=theta2_new;

% pci=pcepoly1d('legendre',epsi(2),20);
pci=pci_epsi2;

% load a.mat
% ahat=a;
a=[]; err=[];
a(1,:)=(gleq(:,3)'*0.5*theta2(:,:))/sum(gleq(:,3)*0.5);
err(1,:)=sum((theta2(:,:)-repmat(a(1,:),size(theta2,1),1)).^2);

for i=2:num_poly
    a(i,:)=((gleq(:,3).*subs(pci(i)))'*theta2(:,:)*0.5)/sum(gleq(:,3).*subs(pci(i)^2)*0.5);
    err(i,:)=sum((theta2(:,:)-repmat(a(1,:),size(theta2,1),1)-subs(pci(2:i))*a(2:i,:)).^2);

% if(max(err(i,:))<1e-8)
%     pci=pci(1:i); 
%     break;
% elseif(max(err(i,:))-max(err(i-1,:))>1e-3)
%     pci=pci(1:i-1);
%     a=a(1:i-1,:);
%     err=err(1:i-1,:);
%     break;
% end
end
% disp(i)

% lim=find(err==max(min(err)));
% [row,~]=find(err==repmat(min(err),20,1));
% a=a(1:min(row),:);
% pci=pci(1:min(row));
% max(max(a-ahat))
% min(min(a))

epsi1=gleq(:,2);
a_new=[];

for i=1:size(a,1)
    a_i=a(i,cump_theta1<=0.9999);
    cump_theta1_i=cump_theta1(cump_theta1<=0.9999);
    a_new(:,i)=interp1(cump_theta1_i',a_i',unifcdf(epsi1,-1,1),'linear','extrap');
end
% 
% for i=1:size(a,1)
% a_new(:,i)=interp1(cump_theta1',a(i,:)',unifcdf(epsi1,-1,1),'linear','extrap');
% end
a=a_new';
% save a.mat a;

% pci_2=pcepoly1d('legendre',epsi(1),20);
pci_2=pci_epsi1;

% load b.mat
% bhat=b;
b=[]; err=[];
b(1,:)=(gleq(:,3)'*0.5*a')/sum(gleq(:,3)*0.5);
err(1,:)=sum((a'-repmat(b(1,:),length(a'),1)).^2);

for i=2:num_poly
    b(i,:)=((gleq(:,3).*subs(pci_2(i)))'*a'*0.5)/sum(gleq(:,3).*subs(pci_2(i)^2)*0.5);
    err(i,:)=sum((a'-repmat(b(1,:),length(a'),1)-subs(pci_2(2:i))*b(2:i,:)).^2);

% if(max(err(i,:))<1e-6)
%     pci_2=pci_2(1:i); 
%     break;
% elseif(max(err(i,:))-max(err(i-1,:))>1)
%     pci_2=pci_2(1:i-1);
%     b=b(1:i-1,:);
%     err=err(1:i-1,:);
%     break;
% end
end
% disp(i)
% lim=find(err==min(err));
% [row,~]=find(err==repmat(min(err),20,1));
% b=b(1:min(row),:);
% min(row)
% max(row)
% pci_2=pci_2(1:min(row));

% save b.mat b;
theta2=pci*(pci_2*b)'; % theta2=b;
% min(min(b))
% max(max(b-bhat))


% cump_new_pci2 = interp1(double(subs(theta2)),unifcdf(epsi2,-1,1),theta2_i','v5cubic','extrap');
cump_new_pci2 = interp1(double(subs(theta2)),unifcdf(epsi2,-1,1),theta2_i','linear','extrap');
figure(31)
drawnow
title('PCE Approximation 2');
plot(theta2_i(1:end-1),diff(cump_new_pci2));
grid on
% plot(theta2_i(1:end-1),diff(cump_new2));
% plot(theta1_i(1:end-1),diff(cump_theta1_i));





theta1_PCE=theta1;
theta2_PCE=theta2;
% size(theta1_PCE)
% size(theta2_PCE)
% int(theta1_PCE,'epsi1',-1,1)
% int(int(theta2_PCE,'epsi1',-1,1),'epsi2',-1,1)

% epsi1=-1+2*rand(1000,1);
% epsi2=-1+2*rand(1000,1);
% figure(2)
% plot(subs(theta1_PCE),subs(theta2_PCE),'*r')

end
