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
beta    = n_eff*k0;
% waist of fundamental mode (vertical)
w0      = 2e-6;
% waist of refractive index profile (vertical)
sig     = 2e-6;
% grating index contrast
dn_g    = 1;            


%% Define target beam (scalar), propagate to grating plane
% grid in grating plane: horizontal cross section at z=0
Lx = 40e-6;     dx = 0.05e-6;
Ly = 40e-6;     dy = 0.05e-6;
x = (-.5*Lx: dx :.5*Lx);
y = (-.5*Ly: dy :.5*Ly);

[xx, yy] = meshgrid(x, y);

%% pump in x
n_in  = [1, 0, 0];
n_out = [0, 0, 1];
E1 = define_target_field(xx, yy, k0, n_out, 'gaussian', 2.5e-6, 50e-6);

figure(1)
pcolor(x, y, abs(E1))    
xlabel('x'), ylabel('y')
shading flat
axis equal
colorbar

% calculate central grating properties (period, rotation, tilt)
[Lam_grat0, alp_grat0, alp_tilt0] = compute_grating_angles(...
    n_out(1), n_out(2), n_out(3), n_in);
Lam_grat0 = Lam_grat0 * lam;
fprintf('grating period = %2.4e\n', Lam_grat0);
fprintf('angle of grating rotation = %2.1f\n', alp_grat0*180/pi);
fprintf('tilt angle = %2.1f\n', alp_tilt0*180/pi);

% local propagation direction of target field, extracted from E-grid
kx = -1i * apply_cen_1st_diff2(E1, dx, 1) ./ E1; kx = real(kx);
ky = -1i * apply_cen_1st_diff2(E1, dy, 2) ./ E1; ky = real(ky);
kz = sqrt(k0^2 - kx.^2 - ky.^2);

[Lam_grat, alp_grat, alp_tilt] = compute_grating_angles(kx, ky, kz, n_in);
Lam_grat = Lam_grat * lam;
 
% =1 horizontal input polarisation; =2 vertical input polarisation
pol = 1;
% transverse prop. constant of pump
ky = 0*Lam_grat;     

[al1, Ex1, Ey1, Ez1] = compute_scattering(...
    Lam_grat, alp_grat, alp_tilt, w0, sig, beta, dn_g, n_eff, ky, pol);

E1norm = sqrt(Ex1.^2 + Ey1.^2 + Ez1.^2);

%% pump in y
n_in  = [0, 1, 0];
n_out = [0, 0, 1];
E2 = define_target_field(xx, yy, k0, n_out, 'gaussian', 2.5e-6, 50e-6);

figure(2)
pcolor(x, y, abs(E1))    
xlabel('x'), ylabel('y')
shading flat
axis equal
colorbar

% calculate central grating properties (period, rotation, tilt)
[Lam_grat0, alp_grat0, alp_tilt0] = compute_grating_angles(...
    n_out(1), n_out(2), n_out(3), n_in);
Lam_grat0 = Lam_grat0 * lam;
fprintf('grating period = %2.4e\n', Lam_grat0);
fprintf('angle of grating rotation = %2.1f\n', alp_grat0*180/pi);
fprintf('tilt angle = %2.1f\n', alp_tilt0*180/pi);

% local propagation direction of target field, extracted from E-grid
kx = -1i * apply_cen_1st_diff2(E2, dx, 1) ./ E2; kx = real(kx);
ky = -1i * apply_cen_1st_diff2(E2, dy, 2) ./ E2; ky = real(ky);
kz = sqrt(k0^2 - kx.^2 - ky.^2);

[Lam_grat, alp_grat, alp_tilt] = compute_grating_angles(kx, ky, kz, n_in);
Lam_grat = Lam_grat * lam;
 
% =1 horizontal input polarisation; =2 vertical input polarisation
pol = 1;
% transverse prop. constant of pump
ky = beta*ones(size(Lam_grat));     

[al2, Ex2, Ey2, Ez2] = compute_scattering(...
    Lam_grat, alp_grat, alp_tilt, w0, sig, beta, dn_g, n_eff, ky, pol);

E2norm = sqrt(Ex2.^2 + Ey2.^2 + Ez2.^2);

%% efficiencies figures
figure(5)
subplot(221)
pcolor(x, y, al1)    
xlabel('x'), ylabel('y')
title('Scattering efficiency al (1/m)')
shading flat
axis equal
colorbar

subplot(223)
pcolor(x, y, al2)    
xlabel('x'), ylabel('y')
title('Scattering efficiency al (1/m)')
shading flat
axis equal
colorbar

figure(6)
subplot(221)
pcolor(x, y, Ex1./E1norm)    
xlabel('x'), ylabel('y')
title('Scattered field pol.: Ex (norm.)')
shading flat
axis equal
colorbar

subplot(222)
pcolor(x, y, Ey1./E1norm)    
xlabel('x'), ylabel('y')
title('Scattered field pol.: Ey (norm.)')
shading flat
axis equal
colorbar

subplot(223)
pcolor(x, y, Ez1./E1norm)    
xlabel('x'), ylabel('y')
title('Scattered field pol.: Ez (norm.)')
shading flat
axis equal
colorbar

figure(7)
subplot(221)
pcolor(x, y, Ex2./E2norm)    
xlabel('x'), ylabel('y')
title('Scattered field pol.: Ex (norm.)')
shading flat
axis equal
colorbar

subplot(222)
pcolor(x, y, Ey2./E2norm)    
xlabel('x'), ylabel('y')
title('Scattered field pol.: Ey (norm.)')
shading flat
axis equal
colorbar

subplot(223)
pcolor(x, y, Ez2./E2norm)    
xlabel('x'), ylabel('y')
title('Scattered field pol.: Ez (norm.)')
shading flat
axis equal
colorbar

%%
% from this, extract dng (to get correct E) - for no pump depletion, this is all
%   NOTE: there would be different possibilities, because the target field is
%   scalar, but the scattered field has polarisation, e.g. trying to match
%   modulus of E-field, match one polarisation, ...
%   NOTE: can also include pump profile here (e.g. Gaussian in transverse
%   direction)

dng1 = abs(E1)./E1norm;
dng1 = dng1/(max(dng1(:)));

figure(9)
subplot(221)
pcolor(x, y, dng1)    
xlabel('x'), ylabel('y')
title('Grating dng (norm.), no pump depletion')
shading flat
axis equal
colorbar

subplot(222)
contour(x, y, dng1)    
xlabel('x'), ylabel('y')
title('Grating dng (norm.), no pump depletion')
shading flat
axis equal
colorbar

dng2 = abs(E2)./E2norm;
dng2 = dng2/(max(dng2(:)));

subplot(223)
pcolor(x, y, dng2)    
xlabel('x'), ylabel('y')
title('Grating dng (norm.), no pump depletion')
shading flat
axis equal
colorbar

subplot(224)
contour(x, y, dng2)    
xlabel('x'), ylabel('y')
title('Grating dng (norm.), no pump depletion')
shading flat
axis equal
colorbar

% to include pump depletion, integrate the grating "loss" and scale dng
% accordingly larger with grating propagation distance to maintain correct
% E field

% loss = cumsum(dng1.^2 .* al1,2) * (x(2)-x(1));
% 
% selpump = 2;
% if (selpump==1)    % case 1: flat pump at input P=1, choose max dng for pump depletion
%     dngbar_max = sqrt(1/max(loss(:,end)));
% 
%     dngbar = dngbar_max*0.99;
%     P = 1 + 0*xx;
%     P = P - dngbar^2*loss;
%     dng_full = real(dng1.*dngbar./sqrt(P));
% end
% if (selpump==2)    % case 2: flat pump at input P=1, choose max dng for fractionial pump depletion
%     dngbar_max = sqrt(1/max(loss(:,end)));
% 
%     dngbar = dngbar_max*sqrt(0.2);
%     P = 1 + 0*xx;
%     P = P - dngbar^2*loss;
%     dng_full = real(dng1.*dngbar./sqrt(P));
% end
% if (selpump==3)    % case 3: optimised pump at input, choose max dng for 1/2 pump depletion
%     dngbar_max = sqrt(1/max(loss(:,end)));
% 
%     dngbar = dngbar_max*0.99;
%     P = (dngbar_max^2*loss(:,end))*ones(1,length(x));
%     P = P - dngbar^2*loss;
%     dng_full = real(dng1.*dngbar./sqrt(P));
% end
% 
% %%
% subplot(223)
% pcolor(x, y, dng_full)   
% xlabel('x'), ylabel('y')
% title('Grating dng, pump depletion')
% shading flat
% axis equal
% colorbar
% 
% subplot(224)
% contour(x, y, dng_full)   
% xlabel('x'), ylabel('y')
% title('Grating dng, pump depletion')
% shading flat
% axis equal
% colorbar
% %%
% figure(8)
% ax1 = axes;
% pcolor(1e6*x, 1e6*y, P);
% xlabel(ax1, 'x/ {\mu}m'), ylabel(ax1, 'y/ {\mu}m')
% shading flat
% a = colorbar;
% ylabel(a,'Intensity/ [a.u.]')
% % ax1.FontSize = 30;
%     
% % combine this grating strength with phase information to produce schematic
% % of the full grating (to compare with the holographic results)
% 
% xxplot = xx;
% yyplot = yy;
% zzplot = 0*xx;
% 
% %alphagrat0
% phi = angle(E2.*exp(-1i*k0*(n_out(1)*xxplot+n_out(2)*yyplot+n_out(3)*zzplot)));
% 
% ng = dng_full.*sin(2*pi/Lam_grat0*(xxplot*cos(alp_grat0)+yyplot*sin(alp_grat0)-zzplot*tan(alp_tilt0)) - phi);
% figure(9)
% subplot(221)
% pcolor(x, y, ng)    
% xlabel('x'), ylabel('y')
% title('Full grating ng')
% shading flat
% axis equal
% colorbar
% 
% x_p = xx(1,:);
% z_p = (-3:0.05:3)*1e-6;
% 
% [xxplot,zzplot] = meshgrid(x_p,z_p);
% yyplot = 0*xxplot;
% 
% phip = phi((length(xx(:,1))+1)/2,:);     % phase of target field at y=0
% phip = ones(length(z_p),1) * phip;
% 
% dng_fullp = dng_full((length(xx(:,1))+1)/2,:);     % dng at y=0
% dng_fullp = ones(length(z_p),1) * dng_fullp;
% 
% ngp = dng_fullp.*sin(2*pi/Lam_grat0*(xxplot*cos(alp_grat0)+yyplot*sin(alp_grat0)-zzplot*tan(alp_tilt0)) - phip).*exp(-(zzplot/sig).^2);
% 
% subplot(222)
% pcolor(x_p,z_p,ngp)    
% xlabel('x'), ylabel('z')
% title('Full grating ng')
% shading flat
% axis equal
% colorbar
% 
% % finally, calculate the scattered field including phase information,
% % propagate to target position, check result
% %%
% f = figure(11);
% % f.Position % this shows the figure position
% ax1 = axes;
% pl1  = pcolor(1e6*x, 1e6*y, ng);
% shading interp
% % clim([-.6, .6])
% 
% ax2  = axes;
% pl2  = pcolor(1e6*x_p, 1e6*z_p, ngp);
% shading interp
% % colormap(ax2, 'jet')
% ylim([-2.5, 2.5])
% % clim([-.6, .6])
% 
% ax1.XTick = [];
% 
% % Define the size of the plot
% ax1_height = .7;
% ax2_height = .12;
% 
% set(ax1, 'Position',[.13 .28 .685 ax1_height]);
% set(ax2, 'Position',[.13 .12 .685 ax2_height]);
% 
% cb1 = colorbar(ax1,'Position',[.83 .28 .04 ax1_height]);
% cb2 = colorbar(ax2,'Position',[.83 .12 .04 ax2_height]);
% 
% xlabel(ax2, 'x/ {\mu}m')
% ylabel(ax1, 'y/ {\mu}m')
% ylabel(cb1, '{\Delta}n_g')
%
% % ax1.FontSize = 30;
% % ax2.FontSize = 30;