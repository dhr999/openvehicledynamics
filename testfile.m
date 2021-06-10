clc
v = functions.vehicle_model;
tic;
input = drive_straight;
[~,vel,alpha] = v.vehicle(input);
toc;

figure(1)
plot(linspace(0,end_time,4000),vel(:,1)*18/5)
hold on
plot(linspace(0,end_time,4000),input(:,1)/20)
grid on
legend('Speed','input')
xlabel('Time')
ylabel('Longitudinal Speed')

figure(2)
plot(linspace(0,end_time,4000),vel(:,2)*18/5)
hold on
plot(linspace(0,end_time,4000),input(:,3))
grid on
legend('Speed','input')
xlabel('Time')
ylabel('Lateral Speed')
% 
% figure(3)
% plot(linspace(0,end_time,4000),alpha(:,1)*180/pi,'LineWidth',2)
% grid on
% xlabel('Time')
% ylabel('Slip angle')

% alpha = linspace(-20*pi/180,20*pi/180,100);
% a = abs(alpha(2) - alpha(1));
% k = linspace(0,1,100);
% T = Tirepacejkacombined_test;
% Fy = zeros(100,100);
% Ca = zeros(100,100);
% 
% for i=1:length(k)
%    [~,fy] = T.tireforce(k(i),alpha,5000);
%    Fy(:,i) = fy;
%    Ca(:,i) = gradient(fy,a);
% end