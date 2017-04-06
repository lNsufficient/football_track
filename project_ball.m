clear;

t0 = 0;
tf = 3;

t = linspace(t0,tf)';

vMax = 61.4; %Fastest ever recorded football shot
v0 = 20; %less than 61,3889
phi = 30*pi/180;
theta = 90*pi/180;
v = v0*[cos(phi)*cos(theta), sin(phi), cos(phi)*sin(theta)];
x0 = 1;
y0 = 0;
z0 = 1000;
s = [x0, y0, z0];

g = 9.82;
Xpos = ones(size(t))*s + t*v +-g*t.^2/2*[0, 1, 0];
Xpos = [Xpos'; ones(1,length(Xpos))];

P = eye(3,4);
x = pflat(P*Xpos);
[cam_cent, princ_axis] = camera_info(P);
princ_axis = princ_axis*100;

figure(1);
clf;
plot3(Xpos(1,:), Xpos(2,:), Xpos(3,:), '.r')
xlabel('X')
ylabel('Y')
zlabel('Z')
hold on;
quiver3(cam_cent(1), cam_cent(2), cam_cent(3), princ_axis(1), princ_axis(2), princ_axis(3))
%quiver3(cam_cent, princ_axis)

figure(2);
clf;
plot(x(1,:), x(2,:), '.r');

