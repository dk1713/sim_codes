% Design of a tilted grating for a 3D arbitrarily shaped output beam with
% controllable polarisation.
%
% Assumptions:
% 1.     target field scalar, from scattering formula match |E| in plane 
% (without the vertical component) to target
% 2.     grating tilt angle is the same everywhere (not correct? local 
% change of grating leads to change of effective tilt angle?)
% 3.     neglect diffraction of pump beam
%
% Design ideas:
% 1.    Conversing two input beams into one square forming a hashtag shape
% grating profile. This might cause some issues with other input pumps
% interact with other respective grating profiles.
% 2.    Having two layers s.t. two input pump doesn't interact with each
% other's grating designs.

%% Parameters
% wavelength
lam     = 780e-9;
k0      = 2*pi/lam;
c       = 3e8;

% effective index of fundamental mode
n_eff   = 1.4635;
% propagation constant of mode
beta    = k0*n_eff;
% waist of output beam
w0      = 3e-6;
% waist of refractive index profile (vertical)
sig     = 3e-6;
% focal distance surface to waist
dist_f  = 50e-6;
% grating index contrast
dn_g    = 1;    

%% Define target beam (scalar), propagate to grating plane
% grid in grating plane: horizontal cross section at z=0
Lx = 40e-6;     dx = 0.05e-6;
Ly = 40e-6;     dy = 0.05e-6;
x = (-.5*Lx: dx :.5*Lx);
y = (-.5*Ly: dy :.5*Ly);

[xx, yy] = meshgrid(x, y);

% direction vector for output beam
n_out = [0, 0, 1];

% Set the pump types with case 1, 2, 3:
%   1. flat pump at input P=1, choose max dng for pump depletion
%   2. flat pump at input P=1, choose max dng for 50% pump depletion (can
%   be adjusted if needed. Change the number in 'sqrt')
%   3. optimised pump at input
setpump = 3;

%% pump in x (pump 1)
disp('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv')
disp('calculating grating profile with pump in x')
n_in  = [1, 0, 0];
% not k0 but k = beta
E1 = define_target_field(xx, yy, beta, n_out, 'gaussian', w0, dist_f);

% calculate central grating properties (period, rotation, tilt)
[Lam_grat0_1, alp_grat0_1, alp_tilt0_1] = compute_grating_angles(...
    n_out(1), n_out(2), n_out(3), n_in);
Lam_grat0_1 = Lam_grat0_1 * lam/n_eff;
fprintf('grating period = %2.2e [m]\n', Lam_grat0_1);
fprintf('angle of grating rotation, theta_rot = %2.1f [deg]\n', alp_grat0_1*180/pi);
fprintf('tilt angle, theta_tilt = %2.1f [deg]\n', alp_tilt0_1*180/pi);

% local propagation direction of target field, extracted from E-grid
kx = -1i * apply_cen_1st_diff2(E1, dx, 1) ./ E1; kx = real(kx);
ky = -1i * apply_cen_1st_diff2(E1, dy, 2) ./ E1; ky = real(ky);
kz = sqrt(beta^2 - kx.^2 - ky.^2);

[Lam_grat, alp_grat, alp_tilt] = compute_grating_angles(kx, ky, kz, n_in);
Lam_grat = Lam_grat * lam/n_eff;
% figure(111);clf;
% contour(Lam_grat,100);
% colorbar
% axis equal

% figure(112);clf;
% contour(180*alp_tilt/pi,100);
% colorbar

% =1 horizontal input polarisation; =2 vertical input polarisation
pol = 1;
% transverse prop. constant of pump
ky = 0*Lam_grat;     

[al1, Ex1, Ey1, Ez1] = compute_scattering(...
    Lam_grat, alp_grat, alp_tilt, w0, sig, beta, dn_g, n_eff, ky, pol);
fprintf('scattering rate = %2.4e\n', norm(al1));
E1norm = sqrt(Ex1.^2 + Ey1.^2 + Ez1.^2);
disp('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

%% efficiency

efficiency    = 100 - 100*exp(-trapz(x, al1));
fprintf('power efficiency    = %2.4e [dB] \n', trapz(x, 10*alpha_ana/log(10)));
fprintf('power scattered out = %2.1f [%%] \n', 100 - 100*10^(-.1*trapz(x, 10*alpha_ana/log(10))));
fprintf('power scattered out = %2.1f [%%] \n', efficiency);

%% pump in y (pump 2)
disp('vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv')
disp('calculating grating profile with pump in y')
n_in  = [0, 1, 0];
E2 = E1;

% calculate central grating properties (period, rotation, tilt)
[Lam_grat0_2, alp_grat0_2, alp_tilt0_2] = compute_grating_angles(...
    n_out(1), n_out(2), n_out(3), n_in);
Lam_grat0_2 = Lam_grat0_2 * lam/n_eff;
fprintf('grating period = %2.2e [m]\n', Lam_grat0_2);
fprintf('angle of grating rotation = %2.1f [deg]\n', alp_grat0_2*180/pi);
fprintf('tilt angle = %2.1f [deg]\n', alp_tilt0_2*180/pi);

% local propagation direction of target field, extracted from E-grid
kx = -1i * apply_cen_1st_diff2(E2, dx, 1) ./ E2; kx = real(kx);
ky = -1i * apply_cen_1st_diff2(E2, dy, 2) ./ E2; ky = real(ky);
kz = sqrt(beta^2 - kx.^2 - ky.^2);

[Lam_grat, alp_grat, alp_tilt] = compute_grating_angles(kx, ky, kz, n_in);
Lam_grat = Lam_grat * lam/n_eff;
 
% pol =1 horizontal input polarisation; =2 vertical input polarisation
pol = 1;
% transverse prop. constant of pump
ky = beta*ones(size(Lam_grat));     

[al2, Ex2, Ey2, Ez2] = compute_scattering(...
    Lam_grat, alp_grat, alp_tilt, w0, sig, beta, dn_g, n_eff, ky, pol);
fprintf('scattering rate = %2.4e\n', norm(al2));
E2norm = sqrt(Ex2.^2 + Ey2.^2 + Ez2.^2);
disp('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

%% efficiencies figures
figure(1); clf;
subplot(221); pcolor(x, y, al1);
xlabel('x'); ylabel('y');
title('Scattering efficiency al for pump 1 (1/m)');
shading flat; axis equal; colorbar;

subplot(223); pcolor(x, y, al2);
xlabel('x'); ylabel('y');
title('Scattering efficiency al for pump 2 (1/m)');
shading flat; axis equal; colorbar;

%% Ex, Ey, Ez components in first grating design (pump in x+ direction)
figure(2); clf;
subplot(221); pcolor(x, y, Ex1./E1norm);
xlabel('x'); ylabel('y');
title('Scattered field pol.: Ex (norm.)');
shading flat; axis equal; colorbar;

subplot(222); pcolor(x, y, Ey1./E1norm);
xlabel('x'); ylabel('y');
title('Scattered field pol.: Ey (norm.)');
shading flat; axis equal; colorbar;

subplot(223); pcolor(x, y, Ez1./E1norm);
xlabel('x'); ylabel('y');
title('Scattered field pol.: Ez (norm.)');
shading flat; axis equal; colorbar;

subplot(224); pcolor(x, y, abs(E1).*Ey1./E1norm);
xlabel('x'); ylabel('y');
title('Scattered field (Ey)');
shading flat; axis equal; colorbar;

%% Ex, Ey, Ez components in first grating design (pump in y+ direction)
figure(3); clf;
subplot(221); pcolor(x, y, Ex2./E2norm)    
xlabel('x'); ylabel('y');
title('Scattered field pol.: Ex (norm.)');
shading flat; axis equal; colorbar;

subplot(222); pcolor(x, y, Ey2./E2norm);
xlabel('x'); ylabel('y');
title('Scattered field pol.: Ey (norm.)');
shading flat; axis equal; colorbar;

subplot(223); pcolor(x, y, Ez2./E2norm);
xlabel('x'); ylabel('y');
title('Scattered field pol.: Ez (norm.)');
shading flat; axis equal; colorbar;

subplot(224); pcolor(x, y, abs(E2).*Ex2./E2norm);
xlabel('x'); ylabel('y');
title('Scattered field (Ex)')
shading flat; axis equal; colorbar;

%% Compute for dn_g for 2 grating designs
% NOTE: there would be different possibilities, because the target field is
% scalar, but the scattered field has polarisation, e.g. trying to match
% modulus of E-field, match one polarisation, ...
% NOTE: can also include pump profile here (e.g. Gaussian in transverse
% direction)

% pump 1 (in x+ axis direction)
dng1 = abs(E1)./E1norm;
dng1 = dng1/(max(dng1(:)));

figure(4); clf;
subplot(221); pcolor(x, y, dng1); 
xlabel('x'); ylabel('y');
title('Grating dng (norm.), no pump(x+) depletion');
shading flat; axis equal; colorbar;

subplot(222); contour(x, y, dng1);
xlabel('x'); ylabel('y');
title('Grating dng (norm.), no pump(x+) depletion');
shading flat; axis equal; colorbar;

% pump 2 (in y+ axis direction)
dng2 = abs(E2)./E2norm;
dng2 = dng2/(max(dng2(:)));

subplot(223); pcolor(x, y, dng2);
xlabel('x'); ylabel('y');
title('Grating dng (norm.), no pump(y+) depletion');
shading flat; axis equal; colorbar;

subplot(224); contour(x, y, dng2);
xlabel('x'); ylabel('y');
title('Grating dng (norm.), no pump(y+) depletion');
shading flat; axis equal; colorbar;

%% Computing for grating designs with different types of pump depletion:
%% For pump in x+ direction. (pump 1)
loss1 = cumsum(dng1.^2 .* al1, 2) * (x(2)-x(1));

if (setpump==1)         % case 1
    dngbar_max = sqrt(1/max(loss1(:,end)));
    dngbar = dngbar_max*0.99;
    P = 1 + 0*xx;
    P = P - dngbar^2*loss1;
    dng_full = real(dng1.*dngbar./sqrt(P));
elseif (setpump==2)     % case 2
    dngbar_max = sqrt(1/max(loss1(:,end)));
    dngbar = dngbar_max*sqrt(0.1);
    P = 1 + 0*xx;
    P = P - dngbar^2*loss1;
    dng_full = real(dng1.*dngbar./sqrt(P));
elseif (setpump==3)     % case 3
    dngbar_max = sqrt(1/max(loss1(:,end)));
    dngbar = dngbar_max*0.99;
    P = (dngbar_max^2*loss1(:,end)).*ones(1,length(x));
    P = P - dngbar^2*loss1;
    dng_full = real(dng1.*dngbar./sqrt(P));
else
    fprintf('Invalid pump type selected for the pump in x');
    return;
end

%% plots for grating profile and pump in x
figure(11); clf;
ax1 = axes;
pcolor(1e6*x, 1e6*y, P);
xlabel(ax1, 'x/ {\mu}m', 'fontsize', 16);
ylabel(ax1, 'y/ {\mu}m', 'fontsize', 16);
shading flat;
axis equal;
a = colorbar;
ylabel(a,'Intensity/ [a.u.]', 'fontsize', 16);
set(gca, 'FontSize', 16);

xxplot = xx; yyplot = yy; zzplot = 0*xx;
% Lam_grat0 = 5.329689101469081e-07;

% alpha_grat0
phi = angle(E1.*exp(-1i*beta/n_eff*(...
    + n_out(1)*xxplot ...
    + n_out(2)*yyplot ...
    + n_out(3)*zzplot)));

% Compute refractive index distribution, ng, in new rotated frame (phi)
ng  = dng_full.*sin(2*pi./Lam_grat0_1.*(...
    xxplot*cos(alp_grat0_1) + yyplot*sin(alp_grat0_1) ...
    - zzplot*tan(alp_tilt0_1)) - phi);

% compute for ng for y=0 line (sideview)
x_p = xx(1,:); z_p = (-3:0.05:3)*1e-6;
[xxplot,zzplot] = meshgrid(x_p,z_p); yyplot = 0*xxplot;

% phase of target field at y=0
phip = phi((length(xx(:,1))+1)/2,:);
phip = ones(length(z_p),1) * phip;

% dng at y=0
dng_fullp = dng_full((length(xx(:,1))+1)/2,:);
dng_fullp = ones(length(z_p),1) * dng_fullp;

ngp = dng_fullp.*sin(2*pi/Lam_grat0_1*(...
    xxplot*cos(alp_grat0_1) + yyplot*sin(alp_grat0_1)...
    - zzplot*tan(alp_tilt0_1)) - phip).*exp(-(zzplot/sig).^2);

% Temp index distribution for superposed case.
ng_temp = ng; ngp_temp = ngp;

% Plot of the top and side view of the 2D planar grating design:
figure(12); clf;
ax1 = axes; pcolor(1e6*x, 1e6*y, ng);
shading interp;
set(gca, 'FontSize', 16);

ax2  = axes; pcolor(1e6*x_p, 1e6*z_p, ngp);
shading interp;
ylim([-3 3]);
set(gca, 'FontSize', 16);

ax1.XTick = []; ax1_height = .68; ax2_height = 6*ax1_height/40;
set(ax1, 'Position',[.21 .26 .50 ax1_height]);
set(ax2, 'Position',[.21 .13 .50 ax2_height]);

cb1 = colorbar(ax1,'Position',[.75 .26 .04 ax1_height]);
colorbar(ax2,'Position',[.75 .13 .04 ax2_height]);

xlabel(ax2, 'x/ {\mu}m', 'fontsize', 16);
ylabel(ax1, 'y/ {\mu}m', 'fontsize', 16);
ylabel(cb1, '{\Delta}n_g', 'fontsize', 16);


% figure(102); clf;
% pcolor(1e6*x, 1e6*y, ng); shading interp; axis equal;
% 
% figure(103); clf;
% pcolor(1e6*x_p, 1e6*z_p, ngp); shading interp; axis equal;

%% cross section
% figure(40); clf;
% plot(1e6*x, ng_temp(401,:))
% 
% figure(41); clf;
% plot(1e6*x, ngp_temp(61,:))

%% cross section overlap
% must run the holo check first so that it has the dataset here.
% ng = ng_temp(401,:);
% ng = ng/max(ng);
% 
% hol1 = squeeze(E1_int(401, :));
% hol1 = hol1/max(hol1);
% 
% figure(50); clf;
% plot( ...
%     1e6*x, ng, ...
%     1e6*x, hol1, '-', ...
%     'linewidth', 2);
% xlabel('x/ {\mu}m', 'fontsize', 16);
% ylabel('normalised {\Delta}n_g', 'fontsize', 16);
% set(gca, 'FontSize', 16);
% xlim([-10, 10]);
% ylim([-1.1 1.1]);
% 
% ngp = ngp_temp(61,:);
% ngp = ngp/max(ngp);
% 
% hol2 = squeeze(E2_int(61, :));
% hol2 = hol2/max(hol2);
% 
% figure(51); clf;
% plot( ...
%     1e6*x, ngp, ...
%     1e6*x, hol2, '-', ...
%     'linewidth', 2);
% xlabel('x/ {\mu}m', 'fontsize', 16);
% ylabel('normalised {\Delta}n_g', 'fontsize', 16);
% set(gca, 'FontSize', 16);
% xlim([-10, 10]);
% ylim([-1.1 1.1]);
%% For pump in y+ direction. (pump 2)
% The pump case set will be same as the pump 1 for preventing confusion.
loss2 = cumsum(dng2.^2 .* al2, 1) * (y(2)-y(1));

if (setpump==1)         % case 1
    dngbar_max = sqrt(1/max(loss2(end,:)));
    dngbar = dngbar_max*0.99;
    P = 1 + 0*yy;
    P = P - dngbar^2*loss2;
    dng_full = real(dng2.*dngbar./sqrt(P));
elseif (setpump==2)     % case 2
    dngbar_max = sqrt(1/max(loss2(end,:)));
    dngbar = dngbar_max*sqrt(0.1);
    P = 1 + 0*yy;
    P = P - dngbar^2*loss2;
    dng_full = real(dng2.*dngbar./sqrt(P));
elseif (setpump==3)     % case 3
    dngbar_max = sqrt(1/max(loss2(end,:)));
    dngbar = dngbar_max*0.99;
    P = (dngbar_max^2*loss2(end,:)).*ones(length(y),1);
    P = P - dngbar^2*loss2;
    dng_full = real(dng2.*dngbar./sqrt(P));
else
    fprintf('Invalid pump type selected for the pump in y');
    return;
end

%% plots for grating profile and pump in y
figure(21); clf;
ax1 = axes;
pcolor(1e6*x, 1e6*y, P);
xlabel(ax1, 'x/ {\mu}m'); ylabel(ax1, 'y/ {\mu}m');
shading flat;
a = colorbar;
ylabel(a,'Intensity/ [a.u.]');
    
xxplot = xx; yyplot = yy; zzplot = 0*xx;

% alpha_grat0
phi = angle(E2.*exp(-1i*beta/n_eff*(...
    + n_out(1)*xxplot ...
    + n_out(2)*yyplot ...
    + n_out(3)*zzplot)));

% Compute refractive index distribution, ng, in new rotated frame (phi)
ng = dng_full.*sin(2*pi/Lam_grat0_2*( ...
    xxplot*cos(alp_grat0_2) + yyplot*sin(alp_grat0_2)...
    - zzplot*tan(alp_tilt0_2)));

% compute for ng for y=0 line (sideview)
x_p = xx(1,:); z_p = (-3:0.05:3)*1e-6;
[xxplot,zzplot] = meshgrid(x_p,z_p); yyplot = 0*xxplot;

% phase of target field at y=0
phip = phi((length(xx(:,1))+1)/2,:);     
phip = ones(length(z_p),1) * phip;

% dng at y=0
dng_fullp = dng_full((length(xx(:,1))+1)/2,:);     
dng_fullp = ones(length(z_p),1) * dng_fullp;

ngp = dng_fullp.*sin(2*pi/Lam_grat0_2*(...
    xxplot*cos(alp_grat0_2) + yyplot*sin(alp_grat0_2)...
    - zzplot*tan(alp_tilt0_2)) - phip).*exp(-(zzplot/sig).^2);

figure(22); clf;
ax1 = axes; pcolor(1e6*x, 1e6*y, ng);
shading interp;
set(gca, 'FontSize', 16);

ax2  = axes; pcolor(1e6*x_p, 1e6*z_p, ngp);
shading interp;
ylim([-3 3]);
set(gca, 'FontSize', 16);

ax1.XTick = []; ax1_height = .68; ax2_height = 6*ax1_height/40;
set(ax1, 'Position',[.21 .26 .50 ax1_height]);
set(ax2, 'Position',[.21 .13 .50 ax2_height]);

cb1 = colorbar(ax1,'Position',[.75 .26 .04 ax1_height]);
colorbar(ax2,'Position',[.75 .13 .04 ax2_height]);

xlabel(ax2, 'x/ {\mu}m', 'fontsize', 16);
ylabel(ax1, 'y/ {\mu}m', 'fontsize', 16);
ylabel(cb1, '{\Delta}n_g', 'fontsize', 16);

%% Plot for superposed grating design
ng_temp     = ng_temp + ng;
ngp_temp    = ngp_temp + ngp;

figure(31); clf;
ax1 = axes; pcolor(1e6*x, 1e6*y, ng_temp);
shading interp;
set(gca, 'FontSize', 16);

ax2 = axes; pcolor(1e6*x_p, 1e6*z_p, ngp_temp);
shading interp;
ylim([-3 3]);
set(gca, 'FontSize', 16);

ax1.XTick = []; ax1_height = .68; ax2_height = 6*ax1_height/40;
set(ax1, 'Position',[.21 .26 .50 ax1_height]);
set(ax2, 'Position',[.21 .13 .50 ax2_height]);

cb1 = colorbar(ax1,'Position',[.75 .26 .04 ax1_height]);
colorbar(ax2,'Position',[.75 .13 .04 ax2_height]);

xlabel(ax2, 'x/ {\mu}m', 'fontsize', 16);
ylabel(ax1, 'y/ {\mu}m', 'fontsize', 16);
ylabel(cb1, '{\Delta}n_g', 'fontsize', 16);