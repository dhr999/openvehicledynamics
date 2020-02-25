clc
v = functions.vehicle_model;
tic;
state = v.vehicle(driver_input);
toc;

figure(1)
plot(linspace(0,end_time,4000),state(:,7)*18/5)
grid on
xlabel('Time')
ylabel('Longitudinal Speed')

figure(2)
plot(linspace(0,end_time,4000),state(:,8)*18/5)
grid on
xlabel('Time')
ylabel('Lateral Speed')