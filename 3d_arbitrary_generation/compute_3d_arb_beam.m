%% Computing grating profile for 3D Arbitrary beam generation
% Code to create tilted Bragg grating profile for 3D arbitrary beam
% generation. 2 steps are needed and this codes is 2nd part out of 2 steps.
%       1. Define the required field shape on the surface or propagate the
%       required field shape onto the surface.
%       2. Propagate the field on the surface onto the grating plane.
%       3. Loop through the vertical layers parallel to the slab mode
%       coming from the initial grating. This considers the model where the
%       dispersion is not assumed.
%
%% Device specification
% Indices [1]
n_air   = 1;
n_clad  = 1.4555;
n_core  = 1.4608;
n_eff   = 1.4635;
dn_g    = 1e-3;

% Heights of the layers [m]
h_clad  = 15e-6;
h_core  = 5e-6;

lam     = 780e-9;
k0      = 2*pi*n_air/lam;
%% Target specification
% distance from the top of the chip.
pos_tar     = [30e-6, 10e-6, 50e-6];
% estimated waist at the surface
waist_tar   = [10e-6, 15e-6];
%% Init
% 1. defining the domain
L_x     = 200e-6;
L_y     = 200e-6;

N_x     = 2^10; 
N_y     = N_x/2^5;

x       = linspace(-.5*L_x, .5*L_x, N_x);
y       = linspace(-.5*L_y, .5*L_y, N_y);

% 2. defining the arbitrary beam
n       = 1;
[xx, yy]    = meshgrid(x, y);
k_cen_x     = k0*pos_tar(1)/sqrt(pos_tar(1)^2+pos_tar(2)^2+pos_tar(3)^2);
k_cen_y     = k0*pos_tar(2)/sqrt(pos_tar(1)^2+pos_tar(2)^2+pos_tar(3)^2);

EE        = exp( -(...
    (.5*(xx - pos_tar(1))./waist_tar(1)).^2     ...
    + (.5*(yy - pos_tar(2))./waist_tar(2)).^2   ...
    ).^n );

figure(1)
pcolor(xx, yy, abs(EE).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal
xline(pos_tar(1), 'r');
yline(pos_tar(2), 'r');

%% Fourier space
dk_x    = 2*pi/L_x;
dk_y    = 2*pi/L_y;
k_x     = (-N_x/2:N_x/2-1) * dk_x;
k_y     = (-N_y/2:N_y/2-1) * dk_y;
[kk_x, kk_y] = meshgrid(k_x, k_y);

% Add in central frequency %** '-/+' here

% computing for phases and 
ky_air  = real( sqrt(k0^2 - k_t.^2) );
ky_clad = real( sqrt(k0_clad^2 - k_t.^2) );
ky_core = real( sqrt(k0_core^2 - k_t.^2) );
phase_air   = asin(k_t/(n_air*k0));
phase_clad  = asin(k_t/(n_clad*k0));
phase_core  = asin(k_t/(n_core*k0));

fprintf('--------------------------------------------------------------\n')
fprintf('Checking for limits of central wavenumber\n')
kk_t    = sqrt((kk_x + k_cen_x).^2 + (kk_y + k_cen_y).^2); %** '+/-' here
kk_z    = real(sqrt(k0^2 - kk_t.^2));
fprintf('k_t(x) range: (%3.3e, %3.3e)',min(k_x+k_cen_x),max(k_x+k_cen_x))
fprintf('\nk_cen_x  = %3.3e\n', k_cen_x)
fprintf('k_t(y) range: (%3.3e, %3.3e)',min(k_y+k_cen_y),max(k_y+k_cen_y))
fprintf('\nk_cen_y  = %3.3e\n', k_cen_y)
fprintf('--------------------------------------------------------------\n')

EE_k    = fftshift(fft2(fftshift(EE))) ...
    .* exp(1i*kk_z.*-1*pos_tar(3)); % propagation back onto the surface.
EE      = fftshift(ifft2(fftshift(EE_k)))...
    .* exp(1i*(k_cen_x.*xx + k_cen_y.*yy)); % add back central frequency.

figure(2)
pcolor(xx, yy, abs(EE).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal
xline(pos_tar(1), 'r');yline(pos_tar(2), 'r');
xline(0, 'k');yline(0, 'k');
xline(L_x/4,'g');yline(L_y/4,'g');xline(-L_x/4,'g');yline(-L_y/4,'g');

%% Select the domain of the grating surface
Lg_x    = 100e-6;
Lg_y    = 100e-6;

Ng_x    = N_x/ceil(L_x/Lg_x);
Ng_y    = N_y/ceil(L_y/Lg_y);

index_x = (0:Ng_x-1) + find(x>-Lg_x/2,1);
index_y = (0:Ng_y-1) + find(y>-Lg_y/2,1);

xg      = x(index_x);
yg      = y(index_y);

[xx_g, yy_g] = meshgrid(xg, yg);
EE_g    = EE(index_y,index_x);

figure(3)
pcolor(xx_g, yy_g, abs(EE_g).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal

%% Propagate the field onto the grating plane
sin_psi_x1 = k_cen_x/k0;
sin_psi_y1 = k_cen_y/k0;

sin_psi_x2 = n_air*sin_psi_x1/n_clad;
sin_psi_y2 = n_air*sin_psi_y1/n_clad;

k_cen_x = k0*n_clad*sin_psi_x2;
k_cen_y = k0*n_clad*sin_psi_y2;

EE_g    = EE_g * exp(-1i*(k_cen_x*xx_g + k_cen_y*yy_g));

EE_g_k  = apply_Fresnel(fftshift(fft2(fftshift(EE_g))),  ...
    .* exp(1i*kk_z.*-1*h_clad); % propagation back onto the surface.
EE      = fftshift(ifft2(fftshift(EE_k)))...
    .* exp(1i*(k_cen_x.*xx + k_cen_y.*yy)); % add back central frequency.

%% Create 2D Bragg grating profile slices along the y-axis