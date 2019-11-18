close all;
clear variables;
clc;
tire = functions.Tirepacejkacombined;
tire.camber =0;
[fx,fy] = tire.tireforce([22 2],78,4000)
% r = 2*pi/180;
% alpha = linspace(-r*2,r*5,10);
% k = linspace(-r*5,r*5,100);
% [fx,fy]=tire.tireforce(k,alpha(1),4000);
% plot(fx/1e3,fy/1e3)
% hold on;grid on;
% xlabel('Fx')
% ylabel('Fy')
% for i = 2:10
%     [fx,fy]=tire.tireforce(k,alpha(i),4000);
%     plot(fx/1e3,fy/1e3)
% end
% [alpha,k] = meshgrid(-r:0.01:r,-r:0.01:r);
% surf(k,alpha,fx)
% ylabel('alpha')
% xlabel('kappa')
% zlabel('Force')
