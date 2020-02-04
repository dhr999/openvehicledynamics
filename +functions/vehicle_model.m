classdef vehicle_model
    methods
        function self = vehicle_model()
            self.m = 1425;
            self.ms = 1050;
            self.a = 1.03;
            self.b = 1.55;
            self.l = 2.5;
            self.tf = 1.54;
            self.tr = 1.52;
            self.Iz = 2500;
            self.Ix = 550;
            self.Iy = 350;
            self.hcg = 0.5;
            self.zcg = 0.1;
            self.Kphi = 76.8e3;
            self.Dphi = 3000e3;
            self.Ktheta = 60e3;
            self.Dtheta = 3000e3;
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
            
            T = [driver_input(1),driver_input(2)];  %Torque to left and right wheel
            
            %% DT-roll pitch model
            function dy = eqs(~,y)
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
        tf;tr %Trackwidth
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
    end
end