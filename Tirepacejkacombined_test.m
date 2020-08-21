classdef Tirepacejkacombined_test  % Based on MF-Tire 6.1 by TNO, The Netherlands
    % Tire Used 205./60R15 91 V
    methods
        function self = Tirepacejkacombined_test()
            self.V0=16.7;
            self.r0=0.3135;
            self.pi0=220000;
            self.Fz0=5500;
            self.cz0=209651;
            self.qreo=0.9974;
            self.Breff=8.386;
            self.Dreff=0.25826;
            self.Freff=0.07394;
            self.pCx1=1.579;
            self.pDx1=1.0422;
            self.pDx2=-00.08285;
            self.pDx3=0;
            self.pEx1=0.11113;
            self.pEx2=0.3143;
            self.pEx3=0;
            self.pEx4=0.001719;
            self.pKx1=21.687;
            self.pKx2=13.728;
            self.pKx3=-0.4098;
            self.pVx1=2.0283e-5;
            self.pVx2=1.0568e-4;
            self.pHx1=2.1615e-4;
            self.pHx2=0.0011598;
            self.ppx1=-0.3485;
            self.ppx2=0.37824;
            self.ppx3=-0.09603;
            self.ppx4=0.06518;
            self.rBx1=13.046;
            self.rBx2=9.718;
            self.rBx3=0;
            self.rCx1=0.9995;
            self.rEx1=-0.4403;
            self.rEx2=-0.4663;
            self.rHx1=-9.968e-5;
            self.pCy1=1.338;
            self.pDy1=0.8785;
            self.pDy2=-0.06452;
            self.pDy3=0;
            self.pEy1=-0.8057;
            self.pEy2=-0.6046;
            self.pEy3=0.09854;
            self.pEy4=-6.697;
            self.pEy5=0;
            self.pKy1=15.324;
            self.pKy2=1.715;
            self.pKy3=0.3695;
            self.pKy4=2.0005;
            self.pKy5=0;
            self.pKy6=-0.8987;
            self.pKy7=-0.23303;
            self.pVy1=-0.00661;
            self.pVy2=0.03592;
            self.pVy3=-0.162;
            self.pVy4=-0.4864;
            self.pHy1=-0.001806;
            self.pHy2=0.00352;
            self.ppy1=-0.6255;
            self.ppy2=-0.06523;
            self.ppy3=-0.16666;
            self.ppy4=0.2811;
            self.ppy5=0;
            self.rBy1=10.622;
            self.rBy2=7.82;
            self.rBy3=0.002037;
            self.rBy4=0;
            self.rCy1=1.0587;
            self.rEy1=0.3148;
            self.rEy2=0.004867;
            self.rHy1=0.009472;
            self.rHy2=0.009754;
            self.rVy1=0.05187;
            self.rVy2=4.853e-4;
            self.rVy3=0;
            self.rVy4=94.63;
            self.rVy5=1.8914;
            self.rVy6=23.8;
            self.pFz1 = 0.7098;
            self.qV1 = 7.742e-4;
            self.camber = -0;
            self.Iw = 1.1;
        end
        function[fx,fy] = tireforce(self,k,alpha1,Fz)
            %         function[fx,fy] = tireforce(self,k,alpha1,Fz)
            dpi = (pi-self.pi0)./self.pi0;
%             cz = self.cz0.*(1+self.pFz1.*dpi);
%             r = self.r0.*(self.qreo + self.qV1.*(self.r0.*omega./self.V0));
%             re = r - (self.Fz0./cz).*(self.Freff.*Fz./self.Fz0 + self.Dreff.*atan(self.Breff.*Fz./self.Fz0));
            Fz0dash = self.Fz0;
            dfz = (Fz - Fz0dash)./Fz0dash;
            camber1 = sin(self.camber);
%             alpha1 = zeros([4,1]);
%            
%             for i = 1:4
%                 if i <= 2
%                    Vy = V(2) + dist(1)*psidot;
%                 elseif i > 2
%                    Vy = V(2) - dist(2)*psidot;
%                 end
%                 Vx = V(1) + ((-1)^i)*psidot*trackw/2;
%                 Vcx = (Vx*cos(delta(i))-Vy*sin(delta(i))); Vcy = Vx*sin(delta(i))+Vy*cos(delta(i)); Vsx = Vcx - re.*omega;
%                 alpha1(i) = atan(Vcy/(Vx+1e-4));
%             end
%             
%             k = -Vsx./abs(Vcx);
            %%Longitudinal Force
            %%Pure Slip
            SVx = Fz.*(self.pVx1 + self.pVx2.*dfz);
            Shx = (self.pHx1 + self.pHx2.*dfz);
            Cx = self.pCx1;
            mux = (self.pDx1 + self.pDx3.*dfz).*(1+self.ppx3.*dpi+self.ppx4.*dpi^2).*(1-self.pDx3.*self.camber^2);
            Dx = mux.*Fz;
            Ex = (self.pEx1 + self.pEx2.*dfz + self.pEx3.*dfz.^2).*(1-self.pEx4.*sign(k));
            Kxk = Fz.*(self.pKx1 + self.pKx2.*dfz).*exp(self.pKx3.*dfz).*(1+self.ppx1.*dpi+self.ppx2.*dpi^2);
            Bx = Kxk ./ (Cx.*Dx);
            keq = k + Shx;
            fx0 = Dx.*sin(Cx.*atan(Bx.*keq-Ex.*(Bx.*keq - atan(Bx.*keq))))+SVx;
            %%Combined Slip
            SHxa = self.rHx1;
            Exa = self.rEx1 + self.rEx2.*dfz;
            Cxa = self.rCx1;
            Bxa = (self.rBx1 + self.rBx3.*self.camber^2).*cos(atan(self.rBx2.*k));
            alphas = alpha1 + SHxa; %alpha1 = -Vcy./|Vcx|
            Gxa0 = cos(Cxa.*atan(Bxa.*SHxa-Exa.*(Bxa.*SHxa-atan(Bxa.*SHxa))));
            Gxa = cos(Cxa.*atan(Bxa.*alphas-Exa.*(Bxa.*alphas-atan(Bxa.*alphas))))./Gxa0;
            fx = Gxa.*fx0;
            %%Latreal Forces
            %Pure Slip
            SVyY = Fz.*(self.pVy3 + self.pVy4.*dfz).*camber1;
            SVy = Fz.*(self.pVy1 + self.pVy2.*dfz) + SVyY;
            Kya = self.pKy1.*Fz0dash.*(1+self.ppy1.*dpi).*(1-self.pKy3.*abs(camber1)).*sin(self.pKy4.*atan((Fz./Fz0dash)./((self.pKy2+self.pKy5.*camber1^2).*(1+self.ppy2.*dpi))));
            SHy = (self.pHy1 + self.pHy2.*dfz);
            ay = alpha1 + SHy;
            Cy = self.pCy1;
            muy = (self.pDy1 + self.pDy2.*dfz).*(1+self.ppy3.*dpi+self.ppy4.*dpi^2).*(1-self.pDy3.*camber1^2);
            Dy = muy.*Fz;
            Ey = (self.pEy1 + self.pEy2.*dfz).*(1+self.pEy5.*camber1^2-(self.pEy3 + self.pEy4.*camber1).*sign(ay));
            By = Kya./(Cy.*Dy);
            fy0 = Dx.*sin(Cy.*atan(By.*ay-Ey.*(By.*ay - atan(By.*ay))))+SVy;
            %Combined Slip
            DVyk = muy.*Fz.*(self.rVy1 + self.rVy2.*dfz + self.rVy3.*camber1).*cos(atan(self.rVy4.*alpha1));
            SVyk = DVyk.*sin(self.rVy5.*atan(self.rVy6.*k));
            SHyk = self.rHy1 + self.rHy2.*dfz;
            Eyk = self.rEy1 + self.rEy2.*dfz;
            Cyk = self.rCy1;
            Byk = (self.rBy1 + self.rBy4.*camber1^2).*cos(atan(self.rBy2.*(alpha1-self.rBy3)));
            ks = k + SHyk;
            Gyk0 = cos(Cyk.*atan(Byk.*SHyk-Eyk.*(Byk.*SHyk-atan(Byk.*SHyk))));
            Gyk = cos(Cyk.*atan(Byk.*ks-Eyk.*(Byk.*ks-atan(Byk.*ks))))./Gyk0;
            fy = Gyk.*fy0+SVyk;
        end
    end
    properties
        V0
        r0
        pi0
        Fz0
        cz0
        qreo
        Breff
        Dreff
        Freff
        pCx1
        pDx1
        pDx2
        pDx3
        pEx1
        pEx2
        pEx3
        pEx4
        pKx1
        pKx2
        pKx3
        pVx1
        pVx2
        pHx1
        pHx2
        ppx1
        ppx2
        ppx3
        ppx4
        rBx1
        rBx2
        rBx3
        rCx1
        rEx1
        rEx2
        rHx1
        pCy1
        pDy1
        pDy2
        pDy3
        pEy1
        pEy2
        pEy3
        pEy4
        pEy5
        pKy1
        pKy2
        pKy3
        pKy4
        pKy5
        pKy6
        pKy7
        pVy1
        pVy2
        pVy3
        pVy4
        pHy1
        pHy2
        ppy1
        ppy2
        ppy3
        ppy4
        ppy5
        rBy1
        rBy2
        rBy3
        rBy4
        rCy1
        rEy1
        rEy2
        rHy1
        rHy2
        rVy1
        rVy2
        rVy3
        rVy4
        rVy5
        rVy6
        pFz1
        qV1
        camber
        Iw
    end
end