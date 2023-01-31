clear all; clc;close all;
%para is the argument of all the parameters in the sytem as a vector
%with following elements: 
%para(1)=p; para(2)=\alpha;para(3)=a; para(4)=\gamma
%para(5)=c1; para(6)=b;
%para(7)=c2;para(8)=c3;para(9)=KN;para(10)=r;para(11)=\beta;para(12)=KI
p=0.8;beta=0.5;alpha=0.4;a=0.8;gamma=0.5;b=4;
c1=0.3;c2=0.4;c3=0.1;
KN=124;KI=100;
r=3;
para=[p,alpha,a,gamma,c1,b,c2,c3,KN,r,beta,KI];
%%
%Here y=(E S N I) 
tspan = [0 200];% time interval
y0 = [1060 10 20 1]; %initial condition y0=(E(0),S(0),N(0),I(0))
[t,y] = ode45(@(t,y) odefcn(t,y,para), tspan, y0);