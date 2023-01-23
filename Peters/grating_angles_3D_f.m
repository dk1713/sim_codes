function [lamgrat,alphagrat,alphatilt] = grating_angles_3D_f(xtarv)
% input: xtrav = [target_x,target_y,target_z]
%        where x=pump propagation, z=transverse (in plane), y=vertical (out of plane)
% output: grating period (units of lambda in material), grating direction, tilt angle

% based on grating_angles_3D_v1.m


lam = 1;            % lengths normalised to wavelength
k0 = 2*pi/lam;      % k in material, assumed equal to beta

xtar = xtarv(1);   % not different coordinate system: here z=transverse, y=vertical
ytar = xtarv(3);
ztar = xtarv(2);

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
alphagrat = atan(kgraty./kgratx);   % grating direction (angle)

% normal to the grating planes

ngratx = koutx-kinx;
ngraty = kouty-kiny;
ngratz = koutz-kinz;
nn = sqrt(ngratx.^2+ngraty.^2+ngratz.^2);
ngratx = ngratx./nn;
ngraty = ngraty./nn;
ngratz = ngratz./nn;
alphatilt = pi/2-acos(ngratz);   % grating tilt angle w.r.t. vertical (=0 vertical, =90 deg horizontal planes)

return
