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
pos_tar     = [20e-6, 10e-6, 50e-6];
% estimated waist at the surface
waist_tar   = [5e-6, 10e-6];
%% Init
% 1. defining the domain
L_x     = 200e-6;
L_y     = 200e-6;

N_x     = 2^10; 
N_y     = N_x/(2^5);

x       = linspace(-.5*L_x, .5*L_x, N_x);
y       = linspace(-.5*L_y, .5*L_y, N_y);

% 2. defining the arbitrary beam
n       = 1;
[xx, yy]    = meshgrid(x, y);
psi_x       = acos(pos_tar(1)/pos_tar(3));
psi_y       = acos(pos_tar(2)/pos_tar(3));

EE        = exp( -(...
    (.5*(xx - pos_tar(1))./waist_tar(1)).^2     ...
    + (.5*(yy - pos_tar(2))./waist_tar(2)).^2   ...
    ).^n ) .*exp(1i*k0.*( xx* cos(psi_x) + yy* cos(psi_y) ));

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

k_cen_y     = -k0*cos(psi_y);
k_cen_x     = k0*cos(psi_x);
EE          = EE .* exp(1i*(k_cen_y.*yy + k_cen_x.*xx));

% Addin central frequency
fprintf('--------------------------------------------------------------\n')
fprintf('Checking for limits of central wavenumber\n')
kk_t    = sqrt((kk_x - k_cen_x).^2 + (kk_y - k_cen_y).^2);
kk_z    = real(sqrt(k0^2 - kk_t.^2));
fprintf('min(k_t) = %3.3e\n',   min(min(kk_t)))
fprintf('k_cen_x  = %3.3e\n',   k_cen_x)
fprintf('k_cen_y  = %3.3e\n',   -k_cen_y)
fprintf('max(k_t) = %3.3e\n',   max(max(kk_t)))
fprintf('--------------------------------------------------------------\n')

EE_k    = fftshift(fft2(fftshift(EE))) .* exp(1i*kk_z.*-1*pos_tar(3));
EE      = fftshift(ifft2(fftshift(EE_k))).* exp(-1i*(k_cen_y.*yy + k_cen_x.*xx));

figure(2)
pcolor(xx, yy, abs(EE).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal
xline(pos_tar(1), 'r');
yline(pos_tar(2), 'r');
xline(0, 'k');
yline(0, 'k');