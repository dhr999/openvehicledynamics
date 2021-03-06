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
            self.r0 = 0.3135;%0.2986
            self.Iw = 1.1;
        end
        function[state,vel,alpha] = vehicle(self, driver_input)
            dist = [self.a  self.b];
            g = 9.81;
            Ixx = self.Ix;
            Iyy = self.Iy;
            Izz = self.Iz;
            Iyz = Iyy - Izz;
            Ixy = Ixx - Iyy;
            h = self.hcg - self.zcg;
            tyre = functions.Tirepacejka;
            
            ti = 0;
            tstep = 0.0025;
            tf = ti + tstep;
            end_time = 10;
            z = 1;
            steps = (end_time - ti)/tstep;
            %% Init Cond
            state = zeros([steps 8]);
            vel = zeros([steps 2]);
            alpha = zeros([steps 4]);
            state(z,7:8) = [50 0]*5/18;     %Intial velocity(x,y) in kmph
            vel(z,:) = state(z,7:8);
            V = state(z,7:8);
            omega = V(1)/self.r0*ones([4 1]);
            Fz = self.m/4*ones([4 1]);
            fx = zeros([4 1]);
            
            while tf<end_time
                T = [0;0;driver_input(z,1);driver_input(z,2)];  %Torque to left and right wheel
                delta = [driver_input(z,3);driver_input(z,4);0;0];
                [~,omega] = ode45(@wheeldyna,[ti tf],omega);
                omega = omega(end,:)';
                [fx,fy,alpha1] = tyre.tireforce(V,omega,Fz,delta,dist,state(z,6));%,self.t/2,state(z,3));%);
                Fx = fx(1)*cos(delta(1)) - fy(1)*sin(delta(1)) + fx(2)*cos(delta(2)) - fy(2)*sin(delta(2)) + fx(3) + fx(4);
                Fy = fy(1)*cos(delta(1)) + fx(1)*sin(delta(1)) + fy(2)*cos(delta(2)) + fx(2)*sin(delta(2)) + fy(3) + fy(4);
                Mz = self.a*(fy(1)*cos(delta(1)) + fx(1)*sin(delta(1)) + fy(2)*cos(delta(2)) + fx(2)*sin(delta(2))) - self.b*(fy(3) + fy(4)) + self.t*(-fx(1)*cos(delta(1)) + fy(1)*sin(delta(1)) + fx(2)*cos(delta(2)) - fy(2)*sin(delta(2)) - fx(3) + fx(4));
                [~,y] = ode45(@bodydyna,[ti tf],(state(z,:))');
                alpha(z,:) = alpha1;
                z=z+1;
                state(z,:) = y(end,:);
                V = state(z,7:8);
                ti=tf;
                tf = ti + tstep;
                vel(z,:) = [state(z,7),state(z,8)]*[cos(state(z,3)),-sin(state(z,3));-sin(state(z,3)),cos(state(z,3))];
            end
            
            %% Wheel Dynamics
            function dy = wheeldyna(~,~)
                dy =  (T - self.r0*fx)/self.Iw;
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
        Iw
    end
end