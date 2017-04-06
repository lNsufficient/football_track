function N = getN(x)
%GETN Summary of this function goes here
%   Detailed explanation goes here

meanx = mean(x(1:2,:),2);

stdx = std(x(1:2,:),0,2);

s = mean(1./stdx);
N = diag([s; s; 1]);
N(1:2,end) = -meanx.*s;
end

