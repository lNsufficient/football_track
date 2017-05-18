clear;
close all;
fig_i = 1;

%% Initialize constants
g = 9.82; %m/s^2
D = [0; -1; 0]; %Direction of gravity

%% Get camera matrix

P = eye(3,4);
K = eye(3);

P_tilde = K\P;

R = P_tilde(1:3,1:3);
t = P_tilde(1:3,4);

%% Get the time difference
delta_t = 1/25; %Framerate of video.
nbr_frames = 9;
timeVec = (0:nbr_frames-1)*delta_t;

%% Generate 3D data

X0 = [-5; 0; 0.01];

V0 = [10; 5; 0];

X = X0*ones(1,length(timeVec)) + V0*timeVec + g/2*D*timeVec.^2;

figure(fig_i);
fig_i = fig_i + 1;
clf;
plot3(X(1,:), X(2,:), X(3,:), '.r')
hold on;
quiver3(0, 0, 0, 0, 0, 1, 5)
xlabel('x');
ylabel('y');
zlabel('z');

X = [X; ones(1,size(X,2))];

%% Project 3D data

x = P*X;


lambda_corr = x(3,:); %The correct values that should be found;

x = pflat(x);
x_noDisturb = x;

maxDev = min(max(x(1,:)) - min(x(1,:)), max(x(2,:)) - min(x(2,:)));
distSize = 3;
disturbations = mvnrnd([0 0],distSize*eye(2), nbr_frames)';
figure(6)
plot(disturbations(1,:), disturbations(2,:),'.')
maxDisturb = max(sqrt(sum(disturbations.^2,1)))
x_perc_disturb = (x_noDisturb-x)./x_noDisturb


useDisturbations = 0;
if ~useDisturbations
    maxDisturb = 1;
end
x(1,:) = x(1,:) + useDisturbations*disturbations(1,:);
x(2,:) = x(2,:) + useDisturbations*disturbations(2,:);





x_tilde = pflat(K\x);
figure(fig_i)
clf;
fig_i = fig_i + 1;
plot(x_tilde(1,:), x_tilde(2,:), '.b');

figure(fig_i)
clf;
fig_i = fig_i + 1;
plot(x(1,:), x(2,:), '.b');
hold on;
plot(x_noDisturb(1,:), x_noDisturb(2,:),'r');



%% Establish b and A matrices. 


%Starting not at x0, but x1:
N = size(X,2) - 1;
x0_tilde = x_tilde(:,1);
% i*delta_t is substituted by t;
%b can now quickly be calculated by:
b = -g/2*R*D*(timeVec(2:end)).^2;
b = b(:);

%but I will also do this inside for loop:
% b_old = b;
% b = b*0;

A = zeros(3*N, 3+1+N);

Rdt = R*delta_t;

for i = 1:N
    curr_start_row = (i-1)*3+1;
    curr_end_row = curr_start_row + 2;
    curr_start_col = i+3+1;
    x_index = i + 1; %because first sample is x0_tilde
    A(curr_start_row:curr_end_row, 1:3) = Rdt*i;
    A(curr_start_row:curr_end_row, 4) = x0_tilde;
    A(curr_start_row:curr_end_row, curr_start_col) = -x_tilde(:,x_index);
end

sol = A\b;
Result_and_correct = [sol, [V0; lambda_corr']]
V0_sol = sol(1:3);
X0_sol = R'*(sol(4)*x0_tilde - t);
Result_and_correctVX = [[X0_sol; V0_sol], [X0; V0]]
G = D*g;
[maxDevLSQ, xp_lsq, Xdist_lsq, X_lsq] = reproj_error(R, t, X0_sol, V0_sol, timeVec, G, x_tilde, X)





%% Minimize reprojection error

%Here it is very important that x_tilde is pflatted, which it should
%already be:
x_tilde = pflat(x_tilde);
%The whole x vector is used, so timeVec should be used.
G = D*g;
cf = @(X, gamma) cfunc(X(1:6), gamma, X(end), x_tilde, R, t, timeVec, G)

gamma_l = 0;
gamma_u = maxDisturb

V0_sol = sol(1:3);
X0_sol = R'*(sol(4)*x0_tilde - t);
s = 0;
lb = -inf*ones(7,1);
lb(end) = -0.1;

X_guess = [X0_sol; V0_sol; s];

TOL = 1e-3;
fun = @(X) X(end); %Trying to minimize s, which will be in position 7 of X.
while ((gamma_u - gamma_l) > TOL)

    gamma = (gamma_l + gamma_u)/2;
    cf_fixGamma = @(X) cf(X, gamma);
    [xsol, fsol,exitflag] = fmincon(fun,X_guess,[],[],[] ,[] ,lb,[],cf_fixGamma);
                           %fmincon(fun,x0,      A, b,Aeq,beq,lb,ub,nonlcon)
    
    cf_test = max(cf_fixGamma(xsol));
    exitflag;
    if (fsol <= 0 && exitflag ~= -2 && exitflag ~= 0)
        X_guess = xsol;
        gamma_u = gamma;
        X_guess(end) = 1;
    else 
        gamma_l = gamma;
    end
end

if X_guess(end) < 1
    disp('Did not change solution')
end
cf_fixS_0 = @(X) cfunc(X(1:6), X(end), 0, x_tilde, R, t, timeVec, G)
X_guess(end) = gamma;
lb(end) = 0;
[X_opt, fsol,exitflag] = fmincon(fun,X_guess,[],[],[] ,[] ,lb,[],cf_fixS_0);

disp('OPTIMAL SOLUTION (REPROJECTION ERROR)')
X_opt = X_guess;
X0_new = X_guess(1:3)
V0_new = X_guess(4:6)

Result_and_correctVX = [[X0; V0], [X0_sol; V0_sol], [X0_new; V0_new]]
[maxDevRep, xp_rep, Xdist_rep, X_rep] = reproj_error(R, t, X0_new, V0_new, timeVec, G, x_tilde, X);
disp('First part done')

maxDevs = [maxDevLSQ maxDevRep; Xdist_lsq, Xdist_rep]
latex = [[X0_sol; V0_sol]',  Xdist_lsq, maxDevLSQ, [X0_new; V0_new]', Xdist_rep, maxDevRep]
latex_corr = [X0', V0']
latex_lsq = [[X0_sol; V0_sol]',  Xdist_lsq, maxDevLSQ]
latex_rep = [[X0_new; V0_new]', Xdist_rep, maxDevRep]
%%
figure(fig_i)
clf;
fig_i = fig_i + 1;
plot(x_tilde(1,:), x_tilde(2,:), 'r');
hold on;
plot(xp_lsq(1,:), xp_lsq(2,:),'.r');
plot(xp_rep(1,:), xp_rep(2,:),'.b');

x_noDist_tilde = pflat(K\x_noDisturb);

figure(fig_i)
clf;
fig_i = fig_i + 1;
plot(x_noDist_tilde(1,:), x_noDist_tilde(2,:),'r');
hold on;
plot(xp_lsq(1,:), xp_lsq(2,:),'.r');
plot(xp_rep(1,:), xp_rep(2,:),'.b');

%% USE CONICS
% Find projected conic:

% Setup matrix in (1.4):

% x_p = x_tilde;
% N_conics = size(x_p,2);
% 
% x_vec = (x_p(1,:))';
% y_vec = (x_p(2,:))';
% x2 = x_vec.^2;
% xy = x_vec.*y_vec;
% y2 = y_vec.^2;
% 
% A_conic = [x2, xy, y2, x_vec, y_vec, ones(N_conics,1)];
% 
% [~, ~, V] = svd(A_conic);
% c_vec = V(:,end);
% a = c_vec(1); b = c_vec(2); c = c_vec(3); d = c_vec(4); e = c_vec(5); f = c_vec(6);
% % Eq (1.3)
% C = [a, b/2, d/2; b/2, c, e/2; d/2, e/2, f];
% 
% figure(fig_i)
% fig_i = fig_i + 1;
% plot(diag(x_p'*C*x_p))

%% 
