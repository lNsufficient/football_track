function H = getH(x1, x2, N1, N2)
%GETH constructs the M matrix used for finding H (so x2 ~ H*x1) and then returns the
%corresponding H. HAS NOT BEEN TESTED USING NORMALIZATION MATRIX YET!
%x1 homogenous coordinates in first image;
%x2 homogenous coordinates in second image;
% x = [x1 x2 x3...; y1 y2 y3...; 1 1 1];
%N normalization matrix, N = d + e, where d = diag matrix and d(3,3) =1
% e is [zeros(3,2) [x0; y0; 0]]

if nargin < 3
    N1 = eye(3);
    N2 = N1;
end

x1_tilde = N1*x1;
x2_tilde = N2*x2;

nbr_samples = size(x1,2);
rows = 3*nbr_samples;
cols = 9+nbr_samples;

M = zeros(rows, cols); 
start_col = 10;
for i = 1:nbr_samples
    start_row = (i-1)*3+1;
    vi = x1_tilde(:,i)';
    M(start_row:start_row+2,1:9) = blkdiag(vi, vi, vi);
    ui = x2_tilde(:,i);
    M(start_row:start_row+2,start_col) = ui;
    start_col = start_col + 1;
end

[~, ~, V] = svd(M);
v = V(:,end);
%If ex2, ex4, ex5 is not working, remove the following line, was recently
%added (and not tested). However, this line should be here, and it should
%not cause any troubles.
v = v/norm(v);

H = zeros(3,3);
H(1,:) = v(1:3);
H(2,:) = v(4:6);
H(3,:) = v(7:9);

H = N2\H*N1;


end