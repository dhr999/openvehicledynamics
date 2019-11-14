% close all;
clear all;
clc;
tire = functions.Tirepacejkacombined;
tire.camber =0;
r = 10*pi/180;
% [fx,fy]=tire.tireforce(0,linspace(-r,r,100),4000);
% plot(linspace(-r,r,100),fy)
% grid on
% hold on
% for i = 1:5
%     [fx,fy]=tire.tireforce(k(i),linspace(-r,r,100),4000);
%     plot(linspace(-r,r,100),fy)
%     legend(k(i))
% end
[alpha,k] = meshgrid(-r:0.01:r,-r:0.01:r);
[fx,fy]=tire.tireforce(k,alpha,4000);
surf(alpha,k,fx)