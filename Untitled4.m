T = Tirepacejkacombined_test;


[alpha,k,load] = meshgrid(-15*pi/180:0.01:15*pi/180,-2:0.01:2,300*9.81:50:1000*9.81);

[fx,~] = T.tireforce(k,alpha,load);