classdef vehicle_model
    methods
        function self = vehicle_model()
            self.m = 1760;
            self.ms = 1590;
            self.a = 1.18;
            self.b = 1.77;
            self.l = 2.5;
            self.t = 1.575;
            self.Iz = 2687.1;
            self.Ix = 894.4;
            self.Iy = 350;
            self.hcg = 0.72;
            self.zcg = 0.15;
            self.Kphi = 189e3;
            self.Dphi = 6e3;
            self.Ktheta = 60e3;
            self.Dtheta = 3000e3;
            self.r0 = 0.3135;
        end
        function[state] = vehicleforce(self, driver_input)
            g = 9.81;
            Ixx = self.Ix;
            Iyy = self.Iy;
            Izz = self.Iz;
            Iyz = Iyy - Izz;
            Ixy = Ixx - Iyy;
            h = self.hcg - self.zcg;
            tyre = functions.Tirepacejkacombined;
            
            ti = 0;
            tstep = 0.001;
            tf = ti + tstep;
            z = 1;
            end_time = 10.0;
            [fx,fy] = 
            while tf<end_time
                T = [0;0;driver_input(z,1);driver_input(z,2)];  %Torque to left and right wheel
                [~,omega] = ode45(,[ti tf],init_cond);
            end
            
            %% Wheel Dynamics
            function dy = wheeldyna(~,y)
                dy =  T - self.ro*fx;;
            end
            %% DT-roll pitch model
            function dy = bodydyna(~,y)
                dy = zeros([8 1]);    %[pitch,roll,yaw,pitch rate,roll rate,yaw rate,Vx,Vy]
                dy(1) = y(4);
                dy(2) = y(5);
                dy(3) = y(6);
                dy(4) = (-self.Ktheta*y(1)-self.Dtheta*y(4)+ h*(self.m*g*sin(y(1))*cos(y(2))-Fx*cos(y(1))*cos(y(2)))+y(6)*(y(6)*sin(y(1))*cos(y(1))*(Ixy+Iyz*cos(y(2))^2)-y(5)*(cos(y(1))^2*Ixx+Iyy*sin(y(2))^2*sin(y(1))^2+Izz*cos(y(2))^2*sin(y(1))^2)-y(4)*(sin(y(1))*cos(y(2))*sin(y(2))*Iyz)))/(Iyy*cos(y(2))^2+Izz*sin(y(2))^2);
                dy(5) = (-self.Kphi*y(2)-self.Dphi*y(5)+h*(Fy*cos(y(2))*cos(y(1))+self.m*g*sin(y(2)))+y(6)*Iyz*(y(6)*sin(y(2))*cos(y(2))*cos(y(1))+y(5)*sin(y(1))*sin(y(2))*cos(y(2)))+y(6)*y(4)*(Iyy*cos(y(2))^2+Izz*sin(y(2))^2))/(Ixx*cos(y(1))^2+Iyy*sin(y(1))^2*sin(y(2))^2+Izz*sin(y(1))^2*cos(y(2))^2);
                dy(6) = (Mz - h*(Fx*sin(y(2))+Fy*sin(y(1))*cos(y(2))))/(Ixx*sin(y(1))^2+cos(y(1))^2*(Iyy*sin(y(2))^2+Izz*cos(y(2))^2));
                dy(7) = y(8)*y(6)+ h*(sin(y(1))*cos(y(2))*(y(6)^2+y(4)^2+y(5)^2) - sin(y(2))*dy(6) - 2*cos(y(2))*y(5)*y(6) - cos(y(1))*cos(y(2))*dy(4) + 2*cos(y(1))*sin(y(2))*y(4)*y(5) + sin(y(1))*sin(y(2))*dy(5)) + Fx/self.m;
                dy(8) = -y(7)*y(6) + h*(-sin(y(1))*cos(y(2))*dy(6) - sin(y(2))*y(6)^2 - 2*cos(y(1))*cos(y(2))*y(4)*y(6) + sin(y(1))*sin(y(2))*y(5)*y(6) - sin(y(2))*y(5)^2 + cos(y(2))*dy(5)) + Fy/self.m;
            end
        end
    end
    properties
        m  %Total Mass of the vehicle
        ms %Sprung Mass of the vehicle
        t %Trackwidth
        l %Wheelbase
        a %Distance of the CG from front axle
        b %Distance of the CG from the rear axle
        hcg %Height of the CG
        zcg %Height of the roll axis in line with the CG
        Ix
        Iy
        Iz
        Ktheta
        Dtheta
        Kphi
        Dphi
        r0
    end
end