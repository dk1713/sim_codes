%% Checking for correct beam generation.
% This script checks for compute_realistic_arb_beam.m such that the
% targeted beam remain in the defined domain. There are other methods to
% ensure this to happen but this is the quickest way. Also needed to see
% the images of the beam at the target plane and grating plane.

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

%% Target specification
lam     = 780e-9;

% distance from the top of the glass.
dist    = 1e-3;

%% Specification for environments
% Wavenumbers
k0      = 2*pi*n_air    /   lam;
k0_clad = 2*pi*n_clad   /   lam;
k0_core = 2*pi*n_core   /   lam;

%% Focussed beam feature
% Define grid points
L_x     = 20e-3;

phi     = 0*pi/180;
x_tar   = [dist*tan(phi), -dist];
waist   = 2e-3;
n = 1;

gaussian = @(x, phi, x_tar, w, n) ...
        exp(-1i * k0 * (  sin(phi)*x  )) .* exp(-( (x-x_tar(1))/w ).^( 2*n ));


% Gaussian beam on the surface of the glass
x      = linspace(-.5*L_x, .5*L_x, 2^14)';
E_tar  = gaussian(x, phi, x_tar(1), waist, n);
%%
% checking the shape
figure(12); clf;
plot(x*1e3, abs(E_tar))
xline(x_tar(1)*1e3, '--r')
xlabel('x / [mm]')
ylabel('|E|')
xlim(x_tar(1)*1e3 + [-waist, waist]*1e3*exp(1));
title('Field on target')

%% Propagation onto the grating plane
k       = linspace(-pi, pi, length(x) )'./(x(2) - x(1));
% fprintf('y_R    = %3.3f [mm]\n', mean(y_R)*1e3)
% fprintf('w_0    = %3.3f [mm]\n', mean(w_0)*1e3)
% fprintf('w_1    = %3.3f [mm]\n', mean(w)*1e3)
% fprintf('L_g    = %3.3f [mm]\n', len_g*1e3)

% central frequency to adjust the resolution in fourier space as max(k) >
% wavenumbers, k's.

k_cen   = k0*sin(phi);
E_tar   = E_tar .* exp(1i*k_cen*x);

%%
fprintf('--------------------------------------------------------------\n')
fprintf('Checking for limits of central wavenumber\n')
k_t     = k - k_cen;
fprintf('min(k_t) = %3.3e\n',   min(k_t))
fprintf('k_cen    = %3.3e\n',   - k0*sin(phi))
fprintf('max(k_t) = %3.3e\n',   max(k_t))
fprintf('--------------------------------------------------------------\n')

ky_air  = real( sqrt(k0^2 - k_t.^2) );
ky_clad = real( sqrt(k0_clad^2 - k_t.^2) );
ky_core = real( sqrt(k0_core^2 - k_t.^2) );

%% Propagation onto the target location
phase_air   = asin(k_t/(n_air*k0));
phase_clad  = asin(k_t/(n_clad*k0));
phase_core  = asin(k_t/(n_core*k0));

Ek_tar  = fftshift( fft( fftshift(E_tar) ) );
Ek_grat = exp(-1i*ky_air*-dist).*Ek_tar;
Ek_grat = apply_Fresnel(Ek_grat, phase_air, phase_clad, n_air, n_clad, 's');
Ek_grat = exp(-1i*ky_clad*-h_clad).*Ek_grat;
Ek_grat = apply_Fresnel(Ek_grat, phase_clad, phase_core, n_clad, n_core, 's');
Ek_grat = exp(-1i*ky_core*-h_core/2).*Ek_grat;

E_grat_shifted  = fftshift(  ifft( fftshift(Ek_grat) )  );
E_grat  = E_grat_shifted .* exp(-1i*k_cen*x);

figure(13); clf;
plot(x*1e3, abs(E_grat))
xline(0, '--r');
xlabel('x / [mm]')
ylabel('|E|')
title('Field on grating')