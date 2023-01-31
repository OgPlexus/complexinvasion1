function dydt = odefcn(t,y,para)
%para is the argument of all the parameters in the sytem as a vector
%with following elements: 
%para(1)=p; para(2)=\alpha;para(3)=a; para(4)=\gamma
%para(5)=c1; para(6)=b;
%para(7)=c2;para(8)=c3;para(9)=KN;para(10)=r;para(11)=beta; para(12)=KI
  dydt = zeros(4,1);
  dydt(1) = para(1).*para(11).*y(4)-para(2).*y(1).*(para(3)./(para(3)+y(3)))+para(4).*y(2);%dE/dt
  dydt(2) = para(2).*y(1).*(para(3)./(para(3)+y(3)))-para(5).*y(3).*y(2)-para(4).*y(2);%dS/dt
  dydt(3)=para(6).*y(3).*(1-(y(3)+para(7).*y(2)+para(8).*y(1))./para(9));%dN/dt
  dydt(4)=para(10).*y(4).*(1-y(4)./para(12))-para(11).*y(4);%dI/dt
end