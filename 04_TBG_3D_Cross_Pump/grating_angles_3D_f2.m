function [lamgrat,alphagrat,alphatilt] = grating_angles_3D_f2(xtar,ytar,ztar)
% input: [target_x,target_y,target_z]
%        where x=pump propagation, y=transverse (in plane), z=vertical (out of plane)
% output: grating period (units of lambda in material), grating direction (rad), tilt angle (rad)

% based on grating_angles_3D_v1.m
% v2: change of coordinates of input, allow for vectors or matrices as input


lam = 1;            % lengths normalised to wavelength
k0 = 2*pi/lam;      % k in material, assumed equal to beta

koutx = k0*xtar./sqrt(xtar.^2+ytar.^2+ztar.^2);
kouty = k0*ytar./sqrt(xtar.^2+ytar.^2+ztar.^2);
koutz = k0*ztar./sqrt(xtar.^2+ytar.^2+ztar.^2);

kinx = k0;
kiny = 0;
kinz = 0;

kgratx = kinx-koutx;    % grating K
kgraty = kiny-kouty;

kgratn = sqrt(kgratx.^2+kgraty.^2); % norm of grating K
lamgrat = 2*pi./kgratn;             % grating period
% alphagrat = atan2(kgraty, kgratx);   % grating direction (angle)
alphagrat = atan(kgraty./kgratx);   % grating direction (angle)

ngratx = koutx-kinx;
ngraty = kouty-kiny;
ngratz = koutz-kinz;
nn = sqrt(ngratx.^2+ngraty.^2+ngratz.^2);
ngratx = ngratx./nn;
ngraty = ngraty./nn;
ngratz = ngratz./nn;
alphatilt = pi/2-acos(ngratz);   % grating tilt angle w.r.t. vertical (=0 vertical, =90 deg horizontal planes)

return
