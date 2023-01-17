function [theta_inc, theta_tilt, theta_grat] = grating_angles(varargin)
% input: xtrav = [target_x,target_y,target_z]
%        where x=pump propagation, z=transverse (in plane), y=vertical (out of plane)
% output: grating period (units of lambda in material), grating direction, tilt angle

% based on grating_angles_3D.m

coor_tar    = varargin{1};   % not different coordinate system: here z=transverse, y=vertical
k0          = varargin{2};
dist        = sqrt(coor_tar(1).^2+coor_tar(2).^2+coor_tar(3).^2);

kout_x  = k0*coor_tar(1)./dist;
kout_y  = k0*coor_tar(2)./dist;
kout_z  = k0*coor_tar(3)./dist;

kin_x   = k0;
kin_y   = 0;
kin_z   = 0;

% normal to the grating planes
kg_x = kin_x-kout_x;    % grating K
kg_y = kin_y-kout_y;
kg_z = kin_z-kout_z;

% grating direction (angle)
theta_inc = atan(kg_y./kg_x);

% tilt angle in original frame
theta_grat = atan(abs(kg_x./kg_z));

% grating tilt angle w.r.t. vertical (=0 vertical, =90 [deg] horizontal)
theta_tilt = asin(abs(kg_z)/sqrt(kg_x.^2+kg_y.^2+kg_z.^2));   

return
