% formula of power loss from tilted grating

% v1: original version only, as in Dom's Opt. Express
% v2: extension to input beam propagating at an angle
% v3: further extension to grating angled with respect to propagation
%     direction of pump


lam = 780e-9; 1.55e-6;          % wavelength
c = 3e8;
k0 = 2*pi./lam;

neff = 1.4635;             % effective index of fundamental mode
beta = neff*k0;         % propagation constant of mode
w0 = 3e-6;              % waist of fundamental mode
sigma = 3e-6;           % waist of refractive index profile
dng = 3e-3;             % grating index contrast

% theta = 30 * pi/180;    % tilt angle of grating
% Lambda =  0.8165*lam/neff;  % grating periodicity
% al_grat = -35.2644*pi/180;     % grating direction in (x,z) plane (i.e. direction of grating K)

% use function to calculate the correct grating angles and period for specific target direction

tarv = [1 1 1];     % direction of target (x,y,z)=(forward,vertical,transverse)
[Lambda0,al_grat,theta] = grating_angles_3D_f(tarv);
Lambda = Lambda0*lam/neff;

fprintf('theta_tilt = %2.2f \n', theta*180/pi);
fprintf('theta_inc  = %2.2f \n', al_grat*180/pi);

w = 1/sqrt(1/w0^2+1/sigma^2);   % effective width of mode and grating


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% pump at an angle inside slab %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kz = beta*(-0.3:0.01:0.3);     % transverse k
kx = sqrt(beta^2-kz.^2);        % longitudinal k
al_prop = asin(kz/beta);       % propagation angle


% incident angle relative to grating planes
% ngrat = [-cos(theta), sin(theta), 0];     % normal to grating plane (x,y,z), then rotate in (x,z) plane by al_grat
% nprop = [cos(al_prop), 0, sin(al_prop)];  % prop. direction

ngratx = -cos(theta)*cos(al_grat);   
ngraty = sin(theta);  
ngratz = -cos(theta)*sin(al_grat);   

npropx = cos(al_prop);  npropy = 0*npropx;   npropz = sin(al_prop);
    
al_inc = acos(-ngratx.*npropx);

figure(12)
plot(al_prop*180/pi, al_inc*180/pi, 'linewidth', 2);
xlabel('input prop angle rel. to x-axis [deg]')
ylabel('incident angle rel. grating slice, \theta_{inc} [deg]')


% define polarisation vectors of s and p, for input (pump) and output (Bragg scattered beam):
% define incident plane, i.e. its normal vector (normal to propagation and normal of grating)
% ninc = ngrat x nprop

nincx = ngraty .* npropz - ngratz.*npropy;
nincy = ngratz .* npropx - ngratx.*npropz;
nincz = ngratx .* npropy - ngraty.*npropx;

nn = sqrt(nincx.^2+nincy.^2+nincz.^2);
nincx = nincx./nn;
nincy = nincy./nn;
nincz = nincz./nn;

% input s polarisation (same as normal to incident plane)

ninsx = nincx; 
ninsy = nincy; 
ninsz = nincz;

% input p pol is normal to propagation direction and s-pol vector, ninp ~ nins x nprop

ninpx = ninsy.*npropz - ninsz.*npropy;
ninpy = ninsz.*npropx - ninsx.*npropz;
ninpz = ninsx.*npropy - ninsy.*npropx;

% propagation direction of scattered beam nscat = nprop - 2*ngrat*(nprop.ngrat)

qprojection = npropx.*ngratx+npropy.*ngraty+npropz.*ngratz;

nscatx = npropx - 2*ngratx.*qprojection;
nscaty = npropy - 2*ngraty.*qprojection;
nscatz = npropz - 2*ngratz.*qprojection;

% figure(30)
% plot(al_prop*180/pi,nscatx,al_prop*180/pi,nscaty,al_prop*180/pi,nscatz)
% xlabel('propagation angle')
% ylabel('x,y,z projections')
% legend('x component','y component','z component')
% title('propagation direction of beam reflected by single grating plane')


% output s polarisation: same as input s polarisation

noutsx = nincx; noutsy = nincy; noutsz = nincz;

% p pol is normal to new propagation direction and s-pol vector, noutp ~ nouts x nscat

noutpx = noutsy.*nscatz - noutsz.*nscaty;
noutpy = noutsz.*nscatx - noutsx.*nscatz;
noutpz = noutsx.*nscaty - noutsy.*nscatx;

% figure(31)
% plot(al_prop*180/pi,ninpx,al_prop*180/pi,ninpy,al_prop*180/pi,ninpz, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('x,y,z projections')
% legend('x component','y component','z component')
% title('p-polarised input')
% 
% figure(32)
% plot(al_prop*180/pi,ninsx,al_prop*180/pi,ninsy,al_prop*180/pi,ninsz, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('x,y,z projections')
% legend('x component','y component','z component')
% title('s-polarised input')
% 
% figure(33)
% plot(al_prop*180/pi,noutpx,al_prop*180/pi,noutpy,al_prop*180/pi,noutpz, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('x,y,z projections')
% legend('x component','y component','z component')
% title('p-polarised output')
% 
% figure(34)
% plot(al_prop*180/pi,noutsx,al_prop*180/pi,noutsy,al_prop*180/pi,noutsz, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('x,y,z projections')
% legend('x component','y component','z component')
% title('s-polarised output')


% input vertical (y) polarisation
%ne = [0, 1, 0];     % direction of polarisation
nex = 0;
ney = 1;
nez = 0;

% project of polarisation on s and p

pys = nex.*ninsx + ney.*ninsy + nez.*ninsz;
pyp = nex.*ninpx + ney.*ninpy + nez.*ninpz;

% input horizontal (x-z) polarisation
%ne = [-sin(al_prop), 0, cos(al_prop)];     % direction of polarisation, orthogonal to nprop
nex = -sin(al_prop);
ney = 0;
nez = cos(al_prop);

% project of polarisation on s and p
pxs = nex.*ninsx + ney.*ninsy + nez.*ninsz;
pxp = nex.*ninpx + ney.*ninpy + nez.*ninpz;

figure(13)
plot( ...
    al_prop*180/pi, pxs.^2, ...
    al_prop*180/pi, pxp.^2, ...
    al_prop*180/pi, pys.^2, '--', ...
    al_prop*180/pi, pyp.^2, '--', 'linewidth', 2)
xlabel('input prop angle rel. to x-axis [deg]');
ylabel('normalised power in');
legend( ...
    'horizontal (s-pol)', 'horizontal (p-pol)', ...
    'vertical (s-pol)', 'vertical (p-pol)')

% figure(14)
% plot(al_prop*180/pi,pxs.^2,al_prop*180/pi,pxp.^2, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('normalised power in');
% legend('z-pol s-componet', 'p-pol');
% title('power in horizontally polarised input (s)')

% figure(15)
% plot(al_prop*180/pi,pys.^2,al_prop*180/pi,pyp.^2, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('normalised power in');
% legend('s-pol', 'p-pol');
% title('power in vertically polarised input (p)')

% output propagation direction and angle (from diffraction)

Kgrat = 2*pi./Lambda;   % grating "propagation" constant

kdiffx = kx - Kgrat.*cos(al_grat);      % diffracted beam
kdiffz = kz - Kgrat.*sin(al_grat);
kdiffy = sqrt(beta^2-kdiffx.^2-kdiffz.^2);

phivar3d = real(asin(kdiffy/beta));

% figure(15)
% plot(al_prop*180/pi,phivar3d*180/pi)
% xlabel('propagation angle')
% ylabel('diffracted beam output angle')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scattering coefficients %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpumpnormal = kx*cos(al_grat) + kz*sin(al_grat); % project pump propgation constant onto grating direction (in x,z plane)

alphas3 = sqrt(2*pi)./w0 ./ (2*cos(al_inc).^2).^2 .* (pi./Lambda).^2 ...
    .* (w./sin(2*theta)).^2 .* (dng./neff).^2;
alphas3 = alphas3 .* exp(-2*(w./sin(2*theta)).^2 ...
    .*(kpumpnormal.*cos(theta).^2-pi./Lambda).^2);
alphas3 = alphas3.*sin(phivar3d);

alphap3 = alphas3 .* cos(2*al_inc).^2;

% figure(26)
% plot(...
%     al_prop*180/pi, alphas3,...
%     al_prop*180/pi, alphap3,...
%     al_prop*180/pi, pxp.^2.*alphap3 + pxs.^2.*alphas3,...
%     al_prop*180/pi, pyp.^2.*alphap3 + pys.^2.*alphas3)
% xlabel('propagation angle')
% ylabel('scattering rate')
% legend('s pol','p pol','z pol, avg','y pol, avg')
% title('includes change in kx, incident angle, scatter power direction')

figure(26)
plot(...
    al_prop*180/pi, alphas3,...
    al_prop*180/pi, alphap3, 'linewidth', 2);
xlabel('input prop angle rel. to x-axis [deg]')
ylabel('scattering rate')
legend('s-pol', 'p-pol')
% title('includes change in kx, incident angle, scatter power direction')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scattered electric field %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Es3 = -1 ./ (2*cos(al_inc).^2) .* w./sin(2*theta) .* dng./neff .* pi^1.5 ./Lambda;
Es3 = Es3 .* exp(-(w./(2*sin(2*theta))).^2 .*(2*kpumpnormal.*cos(theta).^2-2*pi./Lambda).^2);

Ep3 = -cos(2*al_inc) .* Es3 ;

% input vertical (y) polarisation

Eyx = Es3.*pys.*noutsx + Ep3.*pyp.*noutpx;
Eyy = Es3.*pys.*noutsy + Ep3.*pyp.*noutpy;
Eyz = Es3.*pys.*noutsz + Ep3.*pyp.*noutpz;

% input horizontal (x-z) polarisation

Exx = Es3.*pxs.*noutsx + Ep3.*pxp.*noutpx;
Exy = Es3.*pxs.*noutsy + Ep3.*pxp.*noutpy;
Exz = Es3.*pxs.*noutsz + Ep3.*pxp.*noutpz;

% figure(41)
% plot(...
%     al_prop*180/pi, Eyx, ...
%     al_prop*180/pi, Eyy, ...
%     al_prop*180/pi, Eyz, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('Output E-field in x,y,z')
% legend('x component','y component','z component')
% title('Output E-field for vertically polarised input (p)')

% figure(42)
% plot(...
%     al_prop*180/pi, Exx, ...
%     al_prop*180/pi, Exy, ...
%     al_prop*180/pi, Exz, 'linewidth', 2)
% xlabel('input prop angle rel. to x-axis [deg]')
% ylabel('Output E-field in x,y,z')
% legend('x component','y component','z component')
% title('Output E-field for horizontally polarised input (s)')

figure(43)
plot(...
    al_prop*180/pi, sqrt(Exx.^2 + Exy.^2 + Exz.^2), ...
    al_prop*180/pi, sqrt(Eyx.^2 + Eyy.^2 + Eyz.^2), '--', 'linewidth', 2)
xlabel('input prop angle rel. to x-axis [deg]')
ylabel('Electric field norm / [Vm^{-1}]')
legend('horizontal (s)', 'vertical (p)')