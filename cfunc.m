function [c, ceq] = cfunc(X_guess, gamma, s, x, R, t, timeVec, G)
%CFUNC Summary of this function goes here
%   Detailed explanation goes here
%timeVec - row vector

X0 = X_guess(1:3);
V0 = X_guess(4:6);

x1 = x(1,:);
x2 = x(2,:);

N = size(x, 2);
oneVec = ones(1,N);
X0_big = X0*oneVec;
X = X0_big + V0*timeVec+G*0.5*timeVec.^2;
%Skapar hela vektorn:
%     xi, i ökar då radnummer ökar:
v1 = diag((x1'*R(3,:)-oneVec'*R(1,:))*X) + x1'*t(3)-oneVec'*t(1);
v2 = diag((x2'*R(3,:)-oneVec'*R(2,:))*X) + x2'*t(3)-oneVec'*t(2);

denom = R(3,:)*X + t(3);



left = sqrt(v1.^2 + v2.^2);

c = left - gamma*denom' - s;


min_denom = min(denom);
max_c = max(c);
% if (s < 0)
%     disp('=============')
%     s
%     disp('=============')
% end
% s
xp = R*X + t*oneVec;
xp = pflat(xp);
x1p = xp(1,:);
x2p = xp(2,:);

actualGamma = sqrt((x1p - x1).^2 + (x2p-x2).^2);
max_actual_gamma = max(actualGamma);
current_gamma = max(left./denom');
ceq = [];
end

