classdef Tirepacejka
    methods
        function self = Tirepacejka()
            self.Fz0 = 3000;
            self.a0 = 0.08;
            self.b = 0.07;
            self.re = 0.3;
            self.Cy = 1.3;
            self.Ey = -1;
            self.Cz = 2.3;
            self.Ez = -2;
            self.Cx = 1.5;
            self.Ex = -1;
            self.muy0 = 1;
            self.mux0 = 1.26;
            self.c1 = 8;
            self.c2 = 1.33;
            self.c3 = 0.25;
            self.c4 = 0.5;
            self.c5 = 1;
            self.c6 = 0.3;
            self.c7 = 100;
            self.c8 = 15;
            self.c9 = 0.3;
            self.c10 = 0;
            self.c11 = 4;
            self.camber = pi./18;
            self.r = 0.3135;
        end
        function[fx,fy,Mz] = tireforce(self,V,omega,Fz,delta,a,psidot)
            %% Lateral Force
            mux = 1.1;muy=1.3;
            Vcx = V(1); Vcy = V(2); Vsx = Vcx - self.r.*omega;
            alpha = delta - ((Vcy+a.*psidot)./Vcx);
            k = -Vsx./abs(Vcx);
            Cfa0 = 30e3;
            Cfa = self.c1.*self.c2.*self.Fz0.*sin(2.*atan(Fz./(self.c2.*self.Fz0)));
            Cfy = self.c5.*Fz;
            Shy = Cfy./Cfa.*self.camber;
            alpha1 = alpha + Shy;
            alphaeq = (Cfa./Cfa0).*(self.Fz0.*self.muy0./(Fz.*muy)).*alpha1;
            Dy = self.muy0.*self.Fz0;
            By = Cfa0./(Dy.*self.Cy);
            fy0 = Dy.*sin(self.Cy.*atan(By.*alphaeq-self.Ey.*(By.*alphaeq - atan(By.*alphaeq))));
            fy = (muy.*Fz./(self.muy0.*self.Fz0)).*fy0;
            %% Longitudinal force
            Cfk = self.c8.*Fz;
            Cfk0 = 40000;
            keq = (self.mux0.*self.Fz0./(mux.*Fz)).*(Cfk./Cfk0).*k;
            Dx = self.mux0.*self.Fz0;
            Bx = Cfk0./(Dx.*self.Cx);
            fx0 = Dx.*sin(self.Cx.*atan(Bx.*keq-self.Ex.*(Bx.*keq - atan(Bx.*keq))));
            fx = (mux.*Fz)./(self.mux0.*self.Fz0).*fx0;
            %% Self aligning Torque
            Cmy = self.c6.*self.b^2./self.re.*Cfk;
            a = self.a0 .* sqrt(Fz./self.Fz0);
            Cma = self.c4.*a.*Cfa;
            Cma0 = 1.2e3;
            Mzr = (Cmy.*self.camber + Cma.*Shy)./(1+self.c7.*alpha.^2);
            Dz = self.c3.*Dy.*self.a0;
            Bz = -Cma0./self.Cz./Dz;
            Mz0 = Dz.*sin(self.Cz.*atan(Bz.*alphaeq - self.Ez.*(Bz.*alphaeq - atan(Bz.*alphaeq))));
            Mz = (muy.*Fz./(self.muy0.*self.Fz0)).*(Cma./Cma0).*(Cfa0./Cfa).*Mz0 + Mzr;
        end
        function[] = plottire(self,Fz,mux,muy)
            k = (-pi/12:0.001:pi/12);
            alpha = (-pi/12:0.001:pi/12);
            [fx,fy,Mz] = self.tireforce(alpha,k,Fz,mux,muy);
            plot(alpha.*180./pi,30.*Mz,'r')
            hold on
            plot(k.*180./pi,fx,'b')
            plot(alpha.*180./pi,fy,'g')
            grid on
            legend(['30Mz';'  Fx';'  Fy']);
            xlabel('Slip Angle(longitudinal./lateral)(deg)');
            %ylabel('Lateral Force(Fy)');
        end
    end
    properties
        Fz0
        a0
        b
        re
        Cy
        Ey
        Cz
        Ez
        Cx
        Ex
        muy0
        mux0
        c1
        c2
        c3
        c4
        c5
        c6
        c7
        c8
        c9
        c10
        c11
        camber
        r
    end
end