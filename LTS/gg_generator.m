function[ax_max,ax_min,ay_max] = gg_generator()
clearc
C = load("Car_data.mat");
T = Tirepacejkacombined_test;
disp('Generating ggV Diagram .....');
disp('Finding maximum velocity');
vx_max = max_vel();
disp('Finding maximum accleration');
[ax_max,vel] = max_acc();
disp('Finding minimum accleration');
ax_min = min_acc();
disp('Finding maximum lateral accleration');
ay_max = max_ay();
disp('ggV Diagram generated');

    function[vx_max] = max_vel()        %To find maximum velocity of the vehicle
        k = 0.101;
        gear_ratio = min(C.Gear_ratio);
        velrpm = (max(C.E_speed)*pi/30/gear_ratio*C.tire_radius); % in m/s
        if f(velrpm) > 0
            vx_max = velrpm;
        else
            veldrag = fzero(@f,[0,velrpm]);
            vx_max = min(velrpm,veldrag);
        end
        function[value] = f(Ux)
            drag = 0.5*C.rho*C.Area*C.Cd*Ux^2;
            downforce = 0.5*C.rho*C.Area*C.Cl*Ux^2;
            ay = 0;
            w = (1-k)*Ux/C.tire_radius;
            w_eng = min(gear_ratio*(w)*30/pi,18000);
            T_wheel = gear_ratio*interp1(C.E_speed,C.E_torque,w_eng,'spline')/2;
            F_wheel = T_wheel/C.tire_radius;
            value = fzero(@g,0);
            function[value] = g(ax)
                %% Tyre Forces
                Fz21 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 + (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz22 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 - (C.Mass*ay/C.trackwidth*C.hg)/2;
                
                [F21x,~] = T.tireforce(k,0,Fz21);
                [F22x,~] = T.tireforce(k,0,Fz22);
                
                F21x = min(F21x,F_wheel);
                F22x = min(F22x,F_wheel);
                a_x = (F21x + F22x - drag)/C.Mass;
                value = a_x - ax;
            end
        end
    end

    function[ax_max,vel] = max_acc()    %To find maximum longitudinal accleration
        vel = (16:4:vx_max);
        ax_max = zeros(1,length(vel));
        for i = 1:length(vel)
            Ux = vel(i);
            options = optimoptions('fmincon','Display','off');
            [~,fval] = fmincon(@gg_acc,0.075,[],[],[],[],0.05,0.2,[],options);
            ax_max(i) = -1*fval;
        end
        function[ax] = gg_acc(k)
            beta = 0;
            Uy = Ux*tan(beta);
            delta = 0;
            r = 0;
            ay = 0;
            gear = interp1(C.veli,[0,1,2,3,4,5,6,7],Ux,'next');
            gear_ratio = interp1(C.Gear,C.Gear_ratio,gear);
            drag = 0.5*C.rho*C.Area*C.Cd*Ux^2;
            downforce = 0.5*C.rho*C.Area*C.Cl*Ux^2;
            vx = Ux + C.trackwidth/2*r;
            vy = Uy + C.a*r;
            Vx = vx*cos(delta) + vy*sin(delta);
            Vy = vy*cos(delta) - vx*sin(delta);
            alpha_11 = atan(-Vy/Vx);
            vx = Ux - C.trackwidth/2*r;
            Vx = vx*cos(delta) + vy*sin(delta);
            Vy = vy*cos(delta) - vx*sin(delta);
            alpha_12 = atan(-Vy/Vx);
            vy = Uy - C.b*r;
            alpha_22 = atan(-vy/vx);
            vx = Ux + C.trackwidth/2*r;
            alpha_21 = atan(-vy/vx);
            
            w = (1-k)*Ux/C.tire_radius;
            w_eng = min(gear_ratio*(w)*30/pi,18000);
            T_wheel = gear_ratio*interp1(C.E_speed,C.E_torque,w_eng,'spline')/2;
            F_wheel = T_wheel/C.tire_radius;
            ax = -1*fzero(@g,[0,20]);
            
            function[value] = g(ax)
                %% Tyre Forces
                Fz11 = (C.Mass*9.81+downforce)/4 -(C.Mass*ax*C.hg/(C.a+C.b))/2 + (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz12 = (C.Mass*9.81+downforce)/4 -(C.Mass*ax*C.hg/(C.a+C.b))/2 - (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz21 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 + (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz22 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 - (C.Mass*ay/C.trackwidth*C.hg)/2;
                
                if Fz11 <0 || Fz22 <0 || Fz12 < 0 || Fz21 <0    % Penalize the solver for exploring non feasible vector space
                    value = 1e6;
                    return
                end
                
                [F21x,~] = T.tireforce(k,alpha_21,Fz21);
                [F22x,~] = T.tireforce(k,alpha_22,Fz22);
                [~,F11y] = T.tireforce(0,alpha_11,Fz11);
                [~,F12y] = T.tireforce(0,alpha_12,Fz12);
                
                F21x = min(F21x,F_wheel);
                F22x = min(F22x,F_wheel);
                a_x = (F21x + F22x - (F11y+F12y)*sin(delta) - drag)/C.Mass;
                value = a_x - ax;
            end
        end
    end

    function[ax_min] = min_acc()        %To find minimum longitudinal accleration
        ax_min = zeros(1,length(vel));
        for i = 1:length(vel)
            Ux = vel(i);
            options = optimoptions('fmincon','Display','off');
            [~,ax_min(i)] = fmincon(@gg_acc,-0.075,[],[],[],[],-0.2,-0.05,[],options);
        end
        function[ax] = gg_acc(k)
            ay = 0;
            dr = C.brake_distribution;
            drag = 0.5*C.rho*C.Area*C.Cd*Ux^2;
            downforce = 0.5*C.rho*C.Area*C.Cl*Ux^2;
            alpha_11 = 0;
            alpha_12 = 0;
            alpha_22 = 0;
            alpha_21 = 0;
            
            ax = fzero(@g,-20);
            function[value] = g(ax)
                %% Tyre Forces
                Fz11 = (C.Mass*9.81+downforce)/4 -(C.Mass*ax*C.hg/(C.a+C.b))/2 + (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz12 = (C.Mass*9.81+downforce)/4 -(C.Mass*ax*C.hg/(C.a+C.b))/2 - (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz21 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 + (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz22 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 - (C.Mass*ay/C.trackwidth*C.hg)/2;
                
                if Fz11 <0 || Fz22 <0 || Fz12 < 0 || Fz21 <0
                    value = 1e6;
                    return
                end
                
                [F21x,~] = T.tireforce(k,alpha_21,Fz21);
                [F22x,~] = T.tireforce(k,alpha_22,Fz22);
                [F11x,~] = T.tireforce(k,alpha_11,Fz11);
                [F12x,~] = T.tireforce(k,alpha_12,Fz12);
                BF = F11x + F12x;
                BR = F21x + F22x;
                
                index = find([abs(BF),abs(BR)]==min(abs(BF),abs(BR)));
                if index == 2
                    Bf = BR*dr/(1-dr);  %Rear Saturation
                    if abs(Bf) > abs(BF)          %Front Saturation on braking
                        Br = BF*(1-dr)/dr;
                        F21x = Br/2;
                        F22x = Br/2;
                        F11x = BF/2;
                        F12x = BF/2;
                    else
                        F12x = Bf/2;
                        F11x = Bf/2;
                        F21x = BR/2;
                        F22x = BR/2;
                    end
                elseif index == 1
                    Br = BF*dr/(1-dr);  %Front Saturation
                    if abs(Br) > abs(BR)          %Rear Saturation on braking
                        Bf = BR*(1-dr)/dr;
                        F21x = BR/2;
                        F22x = BR/2;
                        F11x = Bf/2;
                        F12x = Bf/2;
                    else
                        F12x = BF/2;
                        F11x = BF/2;
                        F21x = Br/2;
                        F22x = Br/2;
                    end
                else
                    Br = BF*(1-dr)/dr;
                    F11x = BF/2;
                    F12x = BF/2;
                    F21x = Br/2;
                    F22x = Br/2;
                end
                a_x = (F11x + F12x + F21x + F22x - drag)/C.Mass;
                value = a_x - ax;
            end
        end
    end

    function[ay_max] = max_ay()         %To find maximum lateral accleration
        acc_count = zeros(length(vel),20);
        ay_max = zeros(length(vel),length(acc_count));
        for i = 1:length(vel)
            acc_count(i,:) = linspace(ax_min(i),ax_max(i),20);
        end
        acc_count(:,1) = ceil(acc_count(:,1));
        acc_count(:,end) = floor(acc_count(:,end));
        dr = C.brake_distribution;
        for i = 1:size(acc_count,1)
            for j = 1:size(acc_count,2)
                target_acc = acc_count(i,j);
                Ux = vel(i);
                [~,fval] = fmincon(@gg_acc,[0,0,0],[],[],[],[],[-0.1,-0.1,-0.1],[deg2rad(4),deg2rad(25),10],[],optimoptions('fmincon','Display','off'));
                ay_max(i,j) = -1*fval;
            end
        end
        
        function[ay_m] = gg_acc(zz)      %Apply fmincon to this function
            beta=zz(1);delta_r=zz(2);r=zz(3);
            Uy = Ux*tan(beta);
            gear = interp1(C.veli,[0,1,2,3,4,5,6,7],Ux,'next');
            gear_ratio = interp1(C.Gear,C.Gear_ratio,gear);
            drag = 0.5*C.rho*C.Area*C.Cd*Ux^2;
            downforce = 0.5*C.rho*C.Area*C.Cl*Ux^2;
            vx = Ux - C.trackwidth*r;
            vy = Uy + C.a*r;
            Vx = vx*cos(delta_r) + vy*sin(delta_r);
            Vy = vy*cos(delta_r) - vx*sin(delta_r);
            alpha_11 = atan(-Vy/Vx);
            vx = Ux + C.trackwidth*r;
            Vx = vx*cos(delta_r) + vy*sin(delta_r);
            Vy = vy*cos(delta_r) - vx*sin(delta_r);
            alpha_12 = atan(-Vy/Vx);
            vy = Uy - C.b*r;
            alpha_22 = atan(-vy/vx);
            vx = Ux - C.trackwidth*r;
            alpha_21 = atan(-vy/vx);
            
            if alpha_21 <0 || alpha_22 < 0
                ay_m = 1e6;
                return
            end
            
            token = 0;
            input = fmincon(@g,[-0.07,-0.5,1],[],[],[],[],[-0.1,-1,0],[0.1,1,inf],[],optimoptions('fmincon','Display','off'));
            token = 1;
            ay_m = -1*g(input);
            
            function[value] = g(xat)             % To converge ax and ay
                k=xat(1);Tb=xat(2);ay=xat(3);
                ax = target_acc;
                
                %% Tyre Forces
                Fz11 = (C.Mass*9.81+downforce)/4 -(C.Mass*ax*C.hg/(C.a+C.b))/2 + (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz12 = (C.Mass*9.81+downforce)/4 -(C.Mass*ax*C.hg/(C.a+C.b))/2 - (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz21 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 + (C.Mass*ay/C.trackwidth*C.hg)/2;
                Fz22 = (C.Mass*9.81+downforce)/4 +(C.Mass*ax*C.hg/(C.a+C.b))/2 - (C.Mass*ay/C.trackwidth*C.hg)/2;
                
                if Fz11 <0 || Fz22 <0 || Fz12 < 0 || Fz21 <0 || k*Tb<0    %Penalize the solver on exploring non feasible vector space
                    value = 1e6;
                    return
                end
                
                %% Powertrain subsysteam
                if Tb > 0
                    [F21x,F21y] = T.tireforce(k,alpha_21+0.0048,Fz21);
                    [F22x,F22y] = T.tireforce(k,alpha_22+0.0048,Fz22);
                    [F11x,F11y] = T.tireforce(0,alpha_11+0.0048,Fz11);
                    [F12x,F12y] = T.tireforce(0,alpha_12+0.0048,Fz12);
                    if Tb ~= 1
                        throttle_factor = tanh(3.6*Tb);
                    elseif Tb == 1
                        throttle_factor = 1;
                    end
                    w = Ux/C.tire_radius;
                    w_eng = min(gear_ratio*(w)*30/pi,18000);
                    T_wheel = throttle_factor*gear_ratio*interp1(C.E_speed,C.E_torque,w_eng,'spline')/2;
                    F_wheel = T_wheel/C.tire_radius;
                    F21x = min([F21x,F_wheel]);
                    F22x = min([F22x,F_wheel]);
                elseif Tb < 0
                    [F21x,F21y] = T.tireforce(k,alpha_21+0.00749,Fz21);
                    [F22x,F22y] = T.tireforce(k,alpha_22+0.00749,Fz22);
                    [F11x,F11y] = T.tireforce(0,alpha_11+0.00749,Fz11);
                    [F12x,F12y] = T.tireforce(0,alpha_12+0.00749,Fz12);
                    BF = F11x + F12x;
                    BR = F21x + F22x;
                    
                    index = find([BF,BR]==min(BF,BR));
                    if index == 2
                        Bf = BR*dr/(1-dr);  %Rear Saturation
                        if Bf > BF          %Front Saturation on braking
                            Br = BF*(1-dr)/dr;
                            F21x = -Tb*Br/2;
                            F22x = -Tb*Br/2;
                            F11x = -Tb*BF/2;
                            F12x = -Tb*BF/2;
                        else
                            F12x = -Tb*Bf/2;
                            F11x = -Tb*Bf/2;
                            F21x = -Tb*BR/2;
                            F22x = -Tb*BR/2;
                        end
                    elseif index == 1
                        Br = BF*dr/(1-dr);  %Front Saturation
                        if Br > BR          %Rear Saturation on braking
                            Bf = BR*(1-dr)/dr;
                            F21x = -Tb*BR/2;
                            F22x = -Tb*BR/2;
                            F11x = -Tb*Bf/2;
                            F12x = -Tb*Bf/2;
                        else
                            F12x = -Tb*BF/2;
                            F11x = -Tb*BF/2;
                            F21x = -Tb*Br/2;
                            F22x = -Tb*Br/2;
                        end
                    else
                        Br = BF*(1-dr)/dr;
                        F11x = -Tb*BF/2;
                        F12x = -Tb*BF/2;
                        F21x = -Tb*Br/2;
                        F22x = -Tb*Br/2;
                    end
                else
                    [F21x,F21y] = T.tireforce(0,alpha_21,Fz21);
                    [F22x,F22y] = T.tireforce(0,alpha_22,Fz22);
                    [F11x,F11y] = T.tireforce(0,alpha_11,Fz11);
                    [F12x,F12y] = T.tireforce(0,alpha_12,Fz12);
                end
                a_x = ((F11x+F12x)*cos(delta_r) + F21x + F22x - (F11y+F12y)*sin(delta_r) - drag)/C.Mass ;
                a_y = (F21y + F22y + (F11y+F12y)*cos(delta_r) + (F11x+F12x)*sin(delta_r)) /C.Mass ;
                value1 = (ax-a_x)^2 + (ay-a_y)^2;
                
                if ~token
                    value = value1;
                else
                    value = a_y;
                end
            end
        end
    end
end