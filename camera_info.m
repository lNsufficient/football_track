function [cam_cent, principal] = camera_info(P)
%CAMERA_INFO returns camera center and principal axis for P

[~, ~, V] = svd(P);
cam_cent = V(:,end);
if cam_cent(end) == 0
    disp('wierd')
end
cam_cent = pflat(cam_cent);
cam_cent = cam_cent(1:3);

principal = P(end, 1:end-1)';


end

