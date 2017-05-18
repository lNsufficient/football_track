function [maxDev, xp, X_maxDist, X] = reproj_error(R, t, X0, V0, timeVec, G, x, X_corr)
%REPROJ_ERROR Summary of this function goes here
%   Detailed explanation goes here

x1 = x(1,:);
x2 = x(2,:);

N = size(x, 2);
oneVec = ones(1,N);
X0_big = X0*oneVec;


X = X0_big + V0*timeVec+G*0.5*timeVec.^2;

xp = R*X + t*oneVec;
xp = pflat(xp);

x1p = xp(1,:);
x2p = xp(2,:);

actualGamma = sqrt((x1p - x1).^2 + (x2p-x2).^2);
maxDev = max(actualGamma);

X_diffs = X - X_corr(1:3,:);
X_dists = sqrt(sum(X_diffs.^2,1));
X_maxDist = max(X_dists);
end

