function c = cfunc(X, gamma, s, x, R, t, timeVec, G)
%CFUNC Summary of this function goes here
%   Detailed explanation goes here
%timeVec - row vector

X0 = X(1:3);
V0 = X(4:6);

x1 = x(1,:);
x2 = x(2,:);

N = size(x, 2);
X0_big = X0*ones(1,N);
X = X0_big + V0*timeVec-G*timeVec.^2/2;
%Skapar hela vektorn:
%     xi, i ökar då radnummer ökar:
v1 = (x1'*R(3,:)-R(1,:))*X + x1'*t(3)-t(1);
v2 = (x2'*R(3,:)-R(2,:))*X + x2'*t(3)-t(2);

denom = R(3,:)*X + t(3);

left = sqrt(v1.^2 + v2.^2);
c = left - gamma*denom' + s;

end

