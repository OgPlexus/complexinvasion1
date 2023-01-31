clear all; close all
p=0.8;beta=0.5;alpha=0.4;a=0.8;gamma=0.5;b=4;
c1=0.3;c2=0.4;c3=0.1;
KN=77.3073196685;KI=100;
r=3;
I=KI*(r-beta)/r;
%I=0;
para=[p*beta*I,alpha,a,gamma,c1,b,c2,0.1,KN];
%%
[E,S,N]=meshgrid(1000:10:1150,0:10:150,0:10:150);
U= para(1).*para(11).*I-para(2).*E.*(para(3)./(para(3)+N))+para(4).*S;
  V = para(2).*E.*(para(3)./(para(3)+N))-para(5).*N.*S-para(4).*S;
  W =para(6).*N.*(1-(N+para(7).*S+para(8).*E)./para(9));
 L=sqrt(U.^2+V.^2+W.^2);
quiver3(N,E,S,W./L,U./L,V./L,1.5,'k');
axis equal tight
%%
close all
quiver(N(:,:,1),E(:,:,1),W(:,:,1)./L(:,:,1),V(:,:,1)./L(:,:,1),0.5,'k');
axis equal tight

%% N=0, I=0
[E,S]=meshgrid(0:1:25,0:1:25);
U=para(1)-para(2).*E+para(4).*S;
V=para(2).*E-para(4).*S;
L=sqrt(U.^2+V.^2);
quiver(E,S,U./L,V./L,0.5,'k');
axis equal tight
%% I=0 N neq 0
[E,S,N]=meshgrid(0:20:210,0:10:210,0:20:210);
  U = para(1)-para(2).*E.*(para(3)./(para(3)+N))+para(4).*S;
  V = para(2).*E.*(para(3)./(para(3)+N))-para(5).*N.*S-para(4).*S;
  W=para(6).*N.*(1-(N+para(7).*S+para(8).*E)./para(9));
L=sqrt(U.^2+V.^2+W.^2);
quiver3(E,S,N,U./L,V./L,W./N,1.5,'k');
axis equal tight
%%
tspan = [0 200];
%InialCon=[200 200 5;20 0 5;200 0 5;0 200 5;100 100 250];
InialCon=100*rand(100,3);
close all
for i=1:length(InialCon)
y0 = InialCon(i,:);
[t,y] = ode45(@(t,y) odefcn(t,y,para), tspan, y0);
YF(i,:)=y(end,:);
 %plot3(y(:,1),y(:,2),y(:,3),'LineWidth',2); hold on
%text(y0(1),y0(2),y0(3),strcat('P_O=(',num2str(y0),')'),'FontSize',14)
drawnow
end
grid on
plot3(0,0,KN,'*r')
xlabel('$E$','interpreter','latex')
ylabel('$S$','interpreter','latex')
set(gca,'FontSize',14);
zlabel('$N$','interpreter','latex')
text(-40, -10, KN+45,'(0 0 K_N)','FontSize',14)
%%
close all;clc;
plot(t,y(:,2),'LineWidth',2)
hold on
plot(t,y(:,1),'LineWidth',2)
xlabel('$t$','interpreter','latex')
set(gca,'FontSize',14);
legend('S(t)','E(t)','Location', 'Best')
ylabel('Population','interpreter','latex')
%%
close all;clc;
plot(t,y(:,2),'LineWidth',2)
hold on
plot(t,y(:,1),'LineWidth',2)
xlabel('$t$','interpreter','latex')
set(gca,'FontSize',14);
legend('S(t)','E(t)','Location', 'Best')
ylabel('Population','interpreter','latex')
%% Periodic Orbit
close all
plot(y(:,3),y(:,1),'LineWidth',2)
hold on;
plot(7.96958172985,1104.53658291,'*r')%N,E
%plot(7.96958172985,13.9418999488,'*')%N,S
xlabel('$N(t)$','interpreter','latex')
set(gca,'FontSize',14);
%ylabel('$S(t)$','interpreter','latex')
ylabel('$E(t)$','interpreter','latex')
legend('\phi_t(X^1(0))','X_{3}^*','Location', 'Best')
%% 3D orbit
close all;
plot3(y(:,3),y(:,1),y(:,2),'LineWidth',2)
grid on
xlabel('$N(t)$','interpreter','latex')
set(gca,'FontSize',14);
ylabel('$E(t)$','interpreter','latex')
zlabel('$S(t)$','interpreter','latex')
 hold on;plot3(7.96958172985,1104.53658291,13.9418999488,'*r')
 %legend('\phi_t(X(0))','X_{3}^*','Location', 'Best')
 %plot3(0.641123866257,540.361652138,173.3067773,'*')
 
%%
%close all;clc;
figure
plot(t(end-1000:end),y(end-1000:end,4),'LineWidth',2);hold on
plot(t(end-1000:end),y(end-1000:end,2),'LineWidth',2)
plot(t(end-1000:end),y(end-1000:end,3),'LineWidth',2)
xlabel('$t$','interpreter','latex')
set(gca,'FontSize',14);
legend('I(t)','S(t)','N(t)','Location', 'Best')
ylabel('Population','interpreter','latex')
%%
clc; %close figure(3); 
figure(3);plot(y(:,3),y(:,2),'LineWidth',2)
hold on;
%plot(1104.53658291,13.9418999488,'*r')%E,S
%plot(7.96958172985,1104.53658291,'*')%N,E
plot(7.96958172985,13.9418999488,'*')%N,S
xlabel('$N(t)$','interpreter','latex')
set(gca,'FontSize',14);
ylabel('$S(t)$','interpreter','latex')
%ylabel('$E(t)$','interpreter','latex')
legend('\phi_t(X^2(0))','X_{3}^*','Location', 'Best')
%%
close all
plot(t,y(:,3),'LineWidth',2)
hold on;
%plot(7.96958172985,1104.53658291,'*r')%N,E
%plot(7.96958172985,13.9418999488,'*')%N,S
xlabel('$t$','interpreter','latex')
set(gca,'FontSize',14);
%ylabel('$S(t)$','interpreter','latex')
ylabel('$N(t)$','interpreter','latex')
%legend('\phi_t(X^1(0))','X_{3}^*','Location', 'Best')
%%
figure
plot(t,y(:,4),'LineWidth',2);hold on
plot(t,y(:,2),'LineWidth',2)
plot(t,y(:,3),'LineWidth',2)
xlabel('$t$','interpreter','latex')
set(gca,'FontSize',14);
legend('I(t)','S(t)','N(t)','Location', 'Best')
ylabel('Population','interpreter','latex')