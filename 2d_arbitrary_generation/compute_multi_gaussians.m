%% Efficiency of an ion trap by beam tracing method
% Refer to the beam tracing method of derivation of the analytic solution.
% By using the mentioned solution, the efficiency is estimated.

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
n       = 1;
%% Target specification
% distance from the top of the chip.
dist    = 8e-3;
% estimated waist at the surface
waist   = 2e-3;
% x position of the target
xs_tar = [100e-6 -100e-6];

%% Filename
model_spec = [                              ...
    '_n_clad_',     n_clad,                 ...
    '_n_core_',     n_core,                 ...
    '_h_clad_',     h_clad,                 ...
    '_h_core_',     h_core,                 ...
    '_lam_',        lam,                    ...
    '_dist_',       dist,                   ...
    '_waist_',      waist,                  ...
    '_spacing_',    xs_tar(2) - xs_tar(1),  ...
    '_num_g_',      length(xs_tar),         ...
    '_max_dng_',    dn_g,                   ...
    ];

%% Specification for environments
% Wavenumbers
k0      = 2*pi*n_air    /   lam;
k0_clad = 2*pi*n_clad   /   lam;
k0_core = 2*pi*n_core   /   lam;

% Define grid points
L_x     = 3e-3;   
x       = linspace(-.5*L_x, .5*L_x, 2^13)';

%% Focussed beam feature
is_focussed = 1;

%% Computing for field at grating plane and target plane
num_of_ite  = length(xs_tar);
mean_phi    = 0;
size_pro    = 500;
E_pro       = zeros(length(x), size_pro);
E_grat      = zeros(size(x));

for ite = 1:num_of_ite
    x_focus = xs_tar(ite); y_focus = dist;
    phi = atan(x_focus/y_focus);
    
    E_surf = gaussian_on_surface(x, lam, n_air, ...
        [x_focus y_focus], waist, n, is_focussed);
    
    % define k
    k       = linspace(-pi, pi, length(x) )'./(x(2) - x(1));
    k_cen   = k0*sin(phi);
    E_surf = E_surf .* exp(1i*k_cen*x);
    
    % Checking for validity of central k (for coarse grid propagation)
    fprintf('----------------------------------------------------------\n')
    fprintf('Checking for limits of central wavenumber\n')
    k_t     = k - k_cen;
    fprintf('min(k_t) = %3.3e\n',   min(k_t))
    fprintf('k_cen    = %3.3e\n',   - k0*sin(phi))
    fprintf('max(k_t) = %3.3e\n',   max(k_t))
    fprintf('----------------------------------------------------------\n')
    
    ky_air  = real( sqrt(k0^2 - k_t.^2) );
    ky_clad = real( sqrt(k0_clad^2 - k_t.^2) );
    ky_core = real( sqrt(k0_core^2 - k_t.^2) );
    phase_air   = asin(k_t/(n_air*k0));
    phase_clad  = asin(k_t/(n_clad*k0));
    phase_core  = asin(k_t/(n_core*k0));

    Ek_surf  = fftshift( fft( fftshift(E_surf) ) );
    Ek_grat = apply_Fresnel(Ek_surf, phase_air, phase_clad, ...
        n_air, n_clad, 's');
    Ek_grat = exp(-1i*ky_clad*-h_clad).*Ek_grat;
    Ek_grat = apply_Fresnel(Ek_grat, phase_clad, phase_core, ...
        n_clad, n_core, 's');
    Ek_grat = exp(-1i*ky_core*-h_core/2).*Ek_grat;

    E_grat_shifted  = fftshift(  ifft( fftshift(Ek_grat) )  );
    E_grat  = E_grat + E_grat_shifted .* exp(-1i*k_cen*x);
    
    % Computing for propagation from surface to the target plane
    y       = linspace(0, dist+.5e-3, size_pro);
    Ek_pro  = exp(-1i*ky_air*y).*(Ek_surf*ones(size(y)));
    E_pro   = E_pro + fftshift(  ifft( fftshift(Ek_pro,1) ),1  ) ...
        .* (exp(-1i*k_cen*x) * ones(size(y)));
    
    mean_phi = mean_phi + phi;
end

mean_phi = mean_phi/num_of_ite;
phi = mean_phi;

% Checking for propagation code
figure(1)
pcolor(1e3*x, 1e3*y, abs(E_pro)')
xlabel('x /[mm]')
ylabel('y /[mm]')
shading flat
colorbar
hold on
xline(xs_tar(1)*1e3, '--r');
yline(dist*1e3, '--r', 'target plane');
xline(0, '--k');
yline(0, '--k', 'on the surface');
hold off
axis equal

%% Required grating profile calculation
% looping over s.t. dn_g is within the boundary of the physical realisable
% limit that depends on pump depletion and limiting to dn_g < 5e-3

phase   = unwrap(angle(E_grat));
dphase  = k0*n_eff + c_diff(x, phase);
period  = 2*pi./abs(dphase);

E_grat  = exp(1i*k0*n_eff*x) .* E_grat;
Pz_amp  = abs(E_grat).^2;
phase   = angle(E_grat);

% init
eta    = .2;
F       = griddedInterpolant(x, Pz_amp, 'spline');
fun     = @(x) F(x);
C       = 1/integral(fun, min(x), max(x));
max_dng = 3e-3;

% parameters in BTA
Lam     = period;
theta   = .5*(  asin(n_clad/n_core*(1/n_clad*sin(phi))) +.5*pi);
w_0     = 2e-6;
sigma   = 2.5e-6;
beta    = 2*pi*n_eff/lam;
K       = 2*pi./Lam;
w_theta = w_0*sigma/sqrt(w_0^2 + sigma^2)/sin(2*theta);

% Finding optinum dng (for loop)
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
    
    if dng_diff > 1
        eta = eta + .05;
    else
        eta = eta + .01*dng_diff;
    end
end
fprintf('optimum dn_g computed, max_dng = %2.5e \n', max_dng);

%% efficiency
kappa       = dn_gs/n_eff;
w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;

alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * (kappa.^2) / w_0 ./ (Lam.^2);
alpha_ana   = alpha_ana .* sin(phi + .5*pi) / (4*cos(theta)^4);
alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);

efficiency  = 100 - 100*exp(-trapz(x, alpha_ana));
fprintf('power efficiency    = %2.4e [dB] \n', trapz(x, 10*alpha_ana/log(10)));
fprintf('power scattered out = %2.1f [%%] \n', 100 - 100*10^(-.1*trapz(x, 10*alpha_ana/log(10))));
fprintf('power scattered out = %2.1f [%%] \n', efficiency);

%% final figures
% two axis
figure(10)
colororder({'b','r'})
yyaxis left
plot(1e3*x, 1e9*period, 'b')
yline(1e9*lam/(1*cos(phi+.5*pi) + n_eff), '--k', 'period for constant k');
xlabel('x / [mm]')
ylabel('grating period, {\Lambda} / [nm]')

yyaxis right
plot(1e3*x, dn_gs, 'r');
xlabel('x / [mm]')
ylabel('index modulation, {\Delta}n_g')

%% Required .mat files
% Save the data generated
filename    =   ['data/tbg_multi_gauss'  model_spec '.mat'];
save(filename, 'n_eff', 'phi', 'x', 'period', 'phase', 'dn_gs', 'efficiency');