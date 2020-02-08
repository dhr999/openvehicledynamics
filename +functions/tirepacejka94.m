classdef tirepacejka94
    methods
        function self = tirepacejka94()
            self.ya = 1000;
            self.a1 = -22.1 ;
            self.a2 = 1011;
            self.a3 = 1078;
            self.a4 = 1.82;
            self.a5 = 0.208;
            self.a6 = 0;
            self.a7 = -0.354;
            self.a8 = 0.707;
            self.a9 = 0.05;
            self.a10 = 0.01;
            self.a11 = 0;
            self.a12 = 0;
            self.a13 = 0;
            self.camber = -2;    %radians
        end
        
        function [fy,Alpha] = lateralForce(self,Fz,alpha)
            ya = self.ya;
            a1 = self.a1;
            a2 = self.a2;
            a3 = self.a3;
            a4 = self.a4;
            a5 = self.a5;
            a6 = self.a6;
            a7 = self.a7;
            a8 = self.a8;
            a9 = self.a9;
            a10 = self.a10;
            a11 = self.a11;
            a12 = self.a12;
            a13 = self.a13;
            camber = self.camber;
            Fz = Fz/1000;
            alpha = alpha*180/pi;
            D = a1*Fz^2 + a2*Fz;
            E = a6*Fz^2 + a7*Fz + a8;
            C = 2-(2/pi)*asin(ya/D);
            BCD = (a3*Fz^2 + a4*Fz)/exp(a5*Fz);
            B = BCD/(C*D);
            Sh = a9*camber;
            Sv = (a10*Fz^2 + a11*Fz)*camber;
            y = D*sin(C*atan(B*alpha-E*(B*alpha - atan(B*alpha))));
            Alpha = alpha + Sh;
            fy = y + Sv;
        end
        
        function[] = plottire(self,Fz)
            alpha = [-pi/12:0.001:pi/12];
            [fy,Alpha] = self.lateralForce(Fz,alpha);
            plot(Alpha,fy)
            grid on
        end
    end
    properties
        ya % Shape factor [-]
        a1 % Load dependency of lateral friction (*1000) [1/kN]
        a2 % Lateral friction level (*1000) [-]
        a3 % Maximum cornering stiffness [N/deg]
        a4 % Load at maximum cornering stiffness [kN]
        a5 % Camber sensitivity of cornering stiffness
        a6 % Load dependency of curvature factor
        a7 % Curvature factor level
        a8 % Camber sensitivity of horizontal shift
        a9 % Load dependency of horizontal shift
        a10 % Horizontal shift level
        a11 % Combined load and camber sensitivity of vertical shift
        a12 % Load dependency of vertical shift
        a13 % Vertical shift level
        camber
    end
end