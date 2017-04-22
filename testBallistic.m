clear;

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
nbr_frames = 25;
t = (0:nbr_frames)*delta_t;

%% Generate 3D data

X0 = [-5; 0; 100];

V0 = [10; 5; 2];

X = X0 + V0*t + g/2*D*t.^2;

figure(1);
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

figure(2)
plot(x(1,:), x(2,:), '.b');

x_tilde = pflat(K\x);
figure(3)
plot(x_tilde(1,:), x_tilde(2,:), '.b');

%% Establish b and A matrices.
%Starting not at x0, but x1:
N = size(X,2) - 1;
x0_tilde = x_tilde(:,1);
% i*delta_t is substituted by t;
%b can now quickly be calculated by:
b = -g/2*R*D*(t(2:end)).^2;
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