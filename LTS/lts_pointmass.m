clearc
load('trackdata.mat')
load('Car_data.mat')

radius = radius_buddha;
sector_length = 0.25;

u_init = 69.92; %(m/s)
steps = length(radius);
%% Create matrices
vel = zeros(1,steps);
gear = zeros(1,steps);
drag = zeros(1,steps);
downforce = zeros(1,steps);
drive_power = zeros(1,steps);
engine_speed = zeros(1,steps);
max_speed = zeros(1,steps);
max_corner_speed = zeros(1,steps);
elapsed_dist = zeros(1,steps);
sector_time = zeros(1,steps);
longi_acc = zeros(1,steps);
lateral_acc = zeros(1,steps);

%% Intialization
vel(1) = u_init;
gear(1) = 7;
elapsed_dist(1) = 0;
% while vel(end) ~= vel(1)
%     if vel(end) == 0
%         vel(1) = u_init;
%     else
%         vel(1) = vel(end);
%     end
for i = 1:steps-1
    drag(i) = 0.5*rho*Area*Cd*vel(i)^2;
    downforce(i) = 0.5*rho*Area*Cl*vel(i)^2;
    normal_force = downforce(i) + Mass*9.81;
    gear_ratio = interp1(Gear,Gear_ratio,gear(i));
    engine_speed(i) = min(vel(i)/tire_radius*gear_ratio*30/pi,18000);
    drive_power(i) = interp1(E_speed,E_power,engine_speed(i),'spline')*0.93;
    propel_force = min((drive_power(i)/vel(i)),sqrt(abs((Mu*normal_force)^2 - (Mass*vel(i)^2/radius(i))^2))/2) - drag(i);
    max_speed(i) = sqrt(vel(i)^2 + 2*sector_length/Mass*propel_force);
    max_corner_speed(i) = power(((Mu*normal_force)^2/((Mass/radius(i))^2 + (drag(i)/vel(i)^2)^2)),0.25);
    vel(i+1) = min(max_corner_speed(i),max_speed(i));
    elapsed_dist(i+1) = elapsed_dist(i) + sector_length;
    downpoint = interp1([1,2,3,4,5,6,7],downshift_point,gear(i));
    if engine_speed(i) == 18000
        gear(i+1) = min(gear(i) + 1,7);
    elseif engine_speed(i) < downpoint
        gear(i+1) = max(gear(i) - 1,1);
    else
        gear(i+1) = gear(i);
    end
end
for a=steps-2:-1:1
    braked_force = drag(a+1) + sqrt((Mu*(downforce(a+1)+Mass*9.81))^2 - (Mass*vel(a+2)^2/radius(a+1))^2);
    if imag(braked_force) ~= 0
        braked_force = 0;
    end
    braked_vel = sqrt(vel(a+2)^2 + 2*sector_length*braked_force/Mass);
    if braked_vel < vel(a+1)
        c_force = drag(a) + sqrt((Mu*(downforce(a)+Mass*9.81))^2 - (Mass*braked_vel^2/radius(a))^2);
        if imag(c_force) ~= 0
            c_force = 0;
        end
        c = sqrt(braked_vel^2 + 2*sector_length*c_force/Mass);
        if c>vel(a)
            vel(a) = braked_vel;
            downforce(a) = 0.5*rho*Area*Cl*vel(a)^2;
            drag(a) = 0.5*rho*Area*Cd*vel(a)^2;
        else
            vel(a) = c;
            vel(a+1) = braked_vel;
            downforce(a) = 0.5*rho*Area*Cl*vel(a)^2;
            downforce(a+1) = 0.5*rho*Area*Cl*vel(a+1)^2;
            drag(a) = 0.5*rho*Area*Cd*vel(a)^2;
            drag(a+1) = 0.5*rho*Area*Cd*vel(a+1)^2;
        end
    end
end
for i = 1:steps-1
    longi_acc(i) = (vel(i+1)^2 - vel(i)^2)/2/sector_length;
    if longi_acc(i) == 0
        longi_acc(i) = 1e-5;
    end
    sector_time(i) = (-vel(i) + sqrt(vel(i)^2 + 2*longi_acc(i)*sector_length))/longi_acc(i);
    lateral_acc(i) = sign_radius_buddha(i)*vel(i)^2/radius(i);
end
lateral_acc(end) = sign(vel(i))*vel(end)^2/radius(end);

elapsed_time = sum(sector_time);
sprintf("Elapsed time = %f",elapsed_time)

plot(elapsed_dist,vel,'b')
hold on
plot(buddha(:,1),buddha(:,2),'r')
grid on