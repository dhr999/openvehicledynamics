T = Tirepacejkacombined_test;
alpha = 0.0574;
k = 0.02;
tic
[~,fy] = T.tireforce(k,alpha,5000);
toc

tic
fy1 = fittedmodel(k,alpha);
toc