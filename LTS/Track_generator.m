x=0;
for i = 1:length(a(:,1))
for j = x+1:floor(x+4*a(i,1))
sign_radius(j) = a(i,2);
end
x = floor(x+4*a(i,1));
end