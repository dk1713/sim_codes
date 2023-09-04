%% Efficiency of an atom trap by beam tracing method
% Refer to the beam tracing method of derivation of the analytic solution.
% By using the mentioned solution, the efficiency is estimated. Make sure
% to use compute_realistic_check.m first to check if you are working in the
% correct domain or the light doesn't go over the domain.

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
dist    = 50e-3;

% grating length
len_g   = 10e-3;

%% Specification for environments
% Wavenumbers
k0      = 2*pi*n_air    /   lam;
k0_clad = 2*pi*n_clad   /   lam;
k0_core = 2*pi*n_core   /   lam;

%% Focussed beam feature
% Define grid points
L_x     = 30e-3;   
Nx      = 2^14; 
dx      = L_x / Nx;
x       = (-Nx/2:Nx/2-1)' * dx;
y_focus = dist;

% rotated frame
phi     = 10*pi/180;
x_focus = dist*tan(phi);
x_rot   = (x - x_focus)*cos(phi) - (0 - y_focus)*sin(phi);
y_rot   = (x - x_focus)*sin(phi) + (0 - y_focus)*cos(phi);

% Waist on the surface of the glass
w       = .5 * len_g / exp(1);

% % beam waist
% w_0     = sqrt(.5*(w^2 + sqrt(w^4 - 4*(dist*lam/pi/n_air)^2)));

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
xline(0, '--r');
xline(-len_g*.5e3, '--k');
xline(len_g*.5e3, '--k');
xlim(1e3*[-len_g, len_g])
xlabel('x / [mm]')
ylabel('|E|')
title('Field on surface')

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
Ek_tar  = exp(-1i*ky_air*dist).* Ek;

E_tar  = fftshift(  ifft( fftshift(Ek_tar) )  ) .* exp(-1i*k_cen*x);

figure(2)
plot(x*1e3, abs(E_tar))
xline(x_focus*1e3, '--r');
xline(x_focus*1e3 - exp(1)*max(w_0)*1e3, '--k');
xline(x_focus*1e3 + exp(1)*max(w_0)*1e3, '--k');
xlim(x_focus*1e3 + [-max(w_0)*4e3 max(w_0)*4e3])
xlabel('x / [mm]')
ylabel('|E|')
title('Field on target')

%% Propagation onto the grating
Ek_grat = apply_Fresnel(Ek, phase_air, phase_clad, n_air, n_clad, 'p');
Ek_grat = exp(-1i*ky_clad*-h_clad).*Ek_grat;
Ek_grat = apply_Fresnel(Ek_grat, phase_clad, phase_core, n_clad, n_core, 'p');
Ek_grat = exp(-1i*ky_core*-h_core/2).*Ek_grat;

E_grat_shifted  = fftshift(  ifft( fftshift(Ek_grat) )  );
E_grat  = E_grat_shifted .* exp(-1i*k_cen*x);

figure(3)
plot(x*1e3, abs(E_grat))
xline(0, '--r');
xlabel('x / [mm]')
ylabel('|E|')
title('Field on grating')

%%
y       = linspace(-.1*dist, 1.1*dist, 500);
Ek_pro  = exp(-1i*ky_air*y).*(Ek*ones(size(y)));
E_pro   = fftshift(  ifft( fftshift(Ek_pro,1) ),1  ) .* (exp(-1i*k_cen*x) * ones(size(y)));

figure(4)
pcolor(x*1e3, y*1e3, abs(E_pro)')
xlabel('x / [mm]')
ylabel('y / [mm]')
shading flat
colorbar
xline(0, '--r');
yline(0, '--r');
xline(1e3*x_focus, '--k');
yline(1e3*y_focus, '--k');
axis equal

%% Required grating profile calculation
phase   = (unwrap(angle(E_grat_shifted)));
dphase  = k0*n_eff - k_cen + c_diff(x, phase);
period  = 2*pi./abs(dphase);

% E_grat  = exp(1i*k0*n_eff*x) .* E_grat;
Pz_amp  = abs(E_grat).^2;

% looping over s.t. dn_g is within the boundary of the physical realisable
% limit that depends on pump depletion and limiting to dn_g < 5e-3

% init
eta    = .2;
F       = griddedInterpolant(x, Pz_amp, 'spline');
fun     = @(x) F(x);
C       = 1/integral(fun, min(x), max(x));
max_dng = 0;

% parameters in BTA
Lam     = period;
theta   = .5*(  asin(n_clad/n_core*(1/n_clad*sin(phi))) +.5*pi);
psi     = acos(lam./Lam - n_eff);
w_0     = 2e-6;
sigma   = 2.5e-6;
beta    = 2*pi*n_eff/lam;
K       = 2*pi./Lam;
w_theta = w_0*sigma/sqrt(w_0^2 + sigma^2)/sin(2*theta);
    
fprintf('computing for optimum dn_g < %2.4e... \n', dn_g);
while 1e3*abs(dn_g - max_dng) > 1e-4
    if eta < 0
        fprintf('eta = %2.2f has to be positive value!\n', eta);
        break
    end
    
    denu = zeros(size(x));
    for i = 1:length(x)
        denu(i) = 1 - eta * C * integral(fun, -100e-6, x(i));
    end

    dng_amp = sqrt(eta * C * fun(x) ./ denu);
    dn_gs   = (lam*n_core/n_eff/w_theta/pi) * sqrt(w_0/sqrt(2*pi))*dng_amp;
    
    max_dng = max(dn_gs); 
    dng_diff = 1e3*(dn_g - max_dng);
    
    if dng_diff < 0
        eta = eta - .05;
    elseif dng_diff > 1
        eta = eta + .05;
    else
        eta = eta + .1*dng_diff;
    end
    fprintf('new max(n_g) = %2.5e \n', max_dng);
end
fprintf('optimum dn_g computed, max_dng = %2.5e \n', max_dng);

%% efficiency
kappa       = dn_gs/n_eff;
w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;

alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * (kappa.^2) / w_0 ./ (Lam.^2);
alpha_ana   = alpha_ana .* sin(psi) / (4*cos(theta)^4);
alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);

% p-pol
% alpha_ana   = alpha_ana*(cos(2*theta).^2);

efficiency    = 100 - 100*exp(-trapz(x, alpha_ana));
fprintf('power efficiency    = %2.4e [dB] \n', trapz(x, 10*alpha_ana/log(10)));
fprintf('power scattered out = %2.1f [%%] \n', 100 - 100*10^(-.1*trapz(x, 10*alpha_ana/log(10))));
fprintf('power scattered out = %2.1f [%%] \n', efficiency);

%% final figures
figure(5)
plot(1e3*x, 1e9*period)
yline(1e9*lam/(1*cos(phi+.5*pi) + n_eff), '--k', 'constant k');
xlabel('x / [mm]')
ylabel('grating period / [nm]')
% xlim([-1.2 1.2])
title('desired grating period')

figure(6)
plot(1e3*x, dn_gs);
xlabel('x / [mm]')
ylabel('grating strength')
% xlim([-1.2 1.2])
title('desired grating strength')

% Required .mat files
save('data_super_real.mat', 'x',  'period', 'dn_gs', 'efficiency');