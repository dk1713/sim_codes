%% Efficiency of an atom trap by beam tracing method
% Refer to the beam tracing method of derivation of the analytic solution.
% By using the mentioned solution, the efficiency is estimated.

%% Device specification
% Indices [1]
n_air   = 1;
n_clad  = 1.4555;
n_core  = 1.4608;
n_eff   = 1.4635;
dn_g    = 5e-3;

% Heights of the layers [m]
h_clad  = 15e-6;
h_core  = 5e-6;

%% Target specification
lam     = 780e-9;

% distance from the top of the glass.
dist    = 5e-3;

% grating length
len_g   = 2e-3;

%% Specification for environments
% Wavenumbers
k0      = 2*pi*n_air    /   lam;
k0_clad = 2*pi*n_clad   /   lam;
k0_core = 2*pi*n_core   /   lam;

%% Focussed beam feature
% Define grid points
L_x     = 1e-3;   
Nx      = 2^13; 
dx      = L_x / Nx;
x       = (-Nx/2:Nx/2-1)' * dx;
y_focus = dist;

% rotated frame
phi     = 0*pi/180;
x_focus = dist*tan(phi);
x_rot   = (x - x_focus)*cos(phi) - (0 - y_focus)*sin(phi);
y_rot   = (x - x_focus)*sin(phi) + (0 - y_focus)*cos(phi);

% Waist on the surface of the glass
w       = .1e-3;

% Rayleigh range(1)
temp_b  = pi*n_air*w^2/lam;
y_R     = .5*( temp_b - sqrt(temp_b^2 - 4*y_focus^2) );

% beam waist(1)
w_0     = real(sqrt( y_R * lam / (pi * n_air) ));
w_rot   = w_0*sqrt(1 + (y_rot./y_R).^2);

% Radius of curvature(y)
R       = y_rot + y_R.^2./y_rot;

% Gouy phase
eta     = atan(y_rot./y_R);
q_inv   = 1./R - 1i*lam./(n_air*pi*w_rot.^2);


% Gaussian beam on the surface of the glass
E      = w_0./w_rot .*exp( -1i*(k0.*y_rot - eta) -1i*k0 * x_rot.^2 .* q_inv/2 );
%%
% checking the shape
figure(1)
plot(x*1e3, abs(E))
xlabel('x / [mm]')
ylabel('|E|')
title('Field on grating plane')

%% Propagation onto the grating plane
k       = linspace(-pi, pi, length(x) )'./(x(2) - x(1));
fprintf('y_R    = %3.3f [mm]\n', mean(y_R)*1e3)
fprintf('w_0    = %3.3f [mm]\n', mean(w_0)*1e3)
fprintf('w_1    = %3.3f [mm]\n', mean(w)*1e3)
fprintf('L_g    = %3.3f [mm]\n', len_g*1e3)

% central frequency to adjust the resolution in fourier space as max(k) >
% wavenumbers, k's.

k_cen   = k0*sin(phi);
E       = E .* exp(1i*k_cen*x);

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

Ek      = fftshift( fft( fftshift(E) ) );
Ek_grat = exp(-1i*ky_core*h_core/2).*Ek;
Ek_grat = apply_Fresnel(Ek_grat, phase_core, phase_clad, n_core, n_clad, 's');
Ek_grat = exp(-1i*ky_clad*h_clad).*Ek_grat;
Ek_grat = apply_Fresnel(Ek_grat, phase_clad, phase_air, n_clad, n_air, 's');

E_surf  = fftshift(  ifft( fftshift(Ek_grat) )  );

figure(2)
plot(x*1e3, abs(E_surf))
xlabel('x / [mm]')
ylabel('|E|')
title('Field on surface')

Ek_tar  = exp(-1i*ky_air*dist).* Ek_grat;

E_shifted  = fftshift(  ifft( fftshift(Ek_grat) )  );
E_tar  = E_shifted .* exp(-1i*k_cen*x);

figure(3)
plot(x*1e3, abs(E_tar))
xlabel('x / [mm]')
ylabel('|E|')
title('Field on target')

%%
y       = linspace(0, 1.1*dist, 500);
Ek_pro  = exp(-1i*ky_air*y).*(Ek_tar*ones(size(y)));
E_pro   = fftshift(  ifft( fftshift(Ek_pro,1) ),1  ) .* (exp(-1i*k_cen*x) * ones(size(y)));

%%
figure(4); clf;
pcolor(x*1e3, y*1e3, abs(E_pro)')
xlabel('x / [mm]')
ylabel('y / [mm]')
shading flat
colorbar
xlim([-.15 .15])
set(gca, 'FontSize', 16);

e_norm = abs(E_pro)';
save('data/prop_data.mat', 'x', 'y', 'e_norm');