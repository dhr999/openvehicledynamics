function[f,Fc,alpha] = cornerforces(u,Xb,Q,delta,Calpha,Rw,A,fy)
    a = A(1);b = A(2);
    B11 = [[0,0,0,0];[-1*Calpha/u,-1*Calpha*a/u,0,0]];
    B12 = [[0,0,0,0];[-1*Calpha/u,-1*Calpha*a/u,0,0]];
    B13 = [[0,0,0,0];[-1*Calpha/u,1*Calpha*b/u,0,0]];
    B14 = [[0,0,0,0];[-1*Calpha/u,1*Calpha*b/u,0,0]];
    B1 = transpose([transpose(B11),transpose(B12),transpose(B13),transpose(B14)]);
    
    W = [Q(1);delta(1);Q(2);delta(2);Q(3);delta(3);Q(4);delta(4)];    
    B21 = [[1/Rw,0];[0,Calpha]];
    B2 = blkdiag(B21,B21,B21,B21);  %%Assumption all tyres are the same
        
    alpha(:) = delta(1:2) - (Xb(1) + Xb(2)*a)/u;
    alpha(3:4) = delta(3:4) - (Xb(1) - Xb(2)*b)/u;
    D1 = zeros(8,1);
    for i=1:4
       D1(2*i) = fy(i) - Calpha*alpha(i); 
    end

    Lw1 = [cos(delta(1)),-sin(delta(1));sin(delta(1)),cos(delta(1))];
    Lw2 = [cos(delta(2)),-sin(delta(2));sin(delta(2)),cos(delta(2))];
    Lw3 = [cos(delta(3)),-sin(delta(3));sin(delta(3)),cos(delta(3))];
    Lw4 = [cos(delta(4)),-sin(delta(4));sin(delta(4)),cos(delta(4))];
    Lw = blkdiag(Lw1,Lw2,Lw3,Lw4);
        
    f = B1*Xb + B2*W + D1;
    Fc = Lw*f;
    alpha = alpha * 180/pi;
end