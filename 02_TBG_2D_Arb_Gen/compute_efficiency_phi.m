%% Efficiency calculation for phis

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
n       = 10;

% distance from the top of the chip.
dist    = 8e-3;
% estimated focussed area
waist   = 30e-6;
L_x     = 10e-3;
% Define grid points 
x = linspace(-.5*L_x, .5*L_x, 2^14)';
k = linspace(-pi, pi, length(x) )'./(x(2) - x(1));

 %% Specification for environments
% Wavenumbers
k0      = 2*pi*n_air    /   lam;
k0_clad = 2*pi*n_clad   /   lam;
k0_core = 2*pi*n_core   /   lam;

%% Focussed beam feature
% Gaussian function handle
gaussian = @(x, phi, x_tar, w, n) ...
    exp(-1i * k0 * (  sin(phi)*x  )) .* exp(-( (x-x_tar(1))/w ).^( 2*n ));

phis            = linspace(-30, 30, 11)*pi/180;
efficiencies    = zeros(size(phis));
for iter = 1:length(phis)
    fprintf('Interation #%i / %i \n', iter, length(phis));
    % tilt angle of the Gaussian
    phi     = phis(iter);
    % position of the target
    x_tar   = [dist*tan(phi), -dist];

    % Gaussian beam
    E_tar   = gaussian(x, phi, x_tar(1), waist, n);

    %% Propagation onto the grating plane
    k_cen   = k0*sin(phi);
    E_tar = E_tar .* exp(1i*k_cen*x);

    %%
    k_t     = k - k_cen;

    ky_air  = real( sqrt(k0^2 - k_t.^2) );
    ky_clad = real( sqrt(k0_clad^2 - k_t.^2) );
    ky_core = real( sqrt(k0_core^2 - k_t.^2) );

    %%
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

    %% Required grating profile calculation
    phase   = (unwrap(angle(E_grat_shifted)));
    dphase  = k0*n_eff - k_cen + c_diff(x, phase);
    period  = 2*pi./abs(dphase);

    E_grat  = exp(1i*k0*n_eff*x) .* E_grat;
    Pz_amp  = abs(E_grat).^2;

    % looping over s.t. dn_g is within the boundary of the physical realisable
    % limit that depends on pump depletion and limiting to dn_g < 5e-3

    % init
    eta     = .01; % 1e-3
%     eta     = .1; % 3e-3
    F       = griddedInterpolant(x, Pz_amp, 'spline');
    fun     = @(x) F(x);
    C       = 1/integral(fun, min(x), max(x));

    % parameters in BTA
    Lam     = period;
    theta   = .5*(  asin(n_clad/n_core*(1/n_clad*sin(phi))) +.5*pi);
    psi     = acos(lam./Lam - n_eff);
    w_0     = 2e-6;
    sigma   = 2.5e-6;
    beta    = 2*pi*n_eff/lam;
    K       = 2*pi./Lam;
    w_theta = w_0*sigma/sqrt(w_0^2 + sigma^2)/sin(2*theta);
    
    dng_diff = 1;
    fprintf('computing for optimum dn_g < %2.4e... \n', dn_g);
    while 1e3*abs(dng_diff) > 1e-4
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
        dng_diff = dn_g - max_dng;

        eta = eta + 5*dng_diff;
        fprintf('showing max dng = %2.4e \n', max(dn_gs));
    end
    %% efficiency
    kappa       = dn_gs/n_eff;
    w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;

    alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * (kappa.^2) / w_0 ./ (Lam.^2);
    alpha_ana   = alpha_ana .* sin(psi) / (4*cos(theta)^4);
    alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);
    
    efficiencies(iter) = real(100 - 100*exp(-trapz(x, alpha_ana)));
end
%%
phi_deg = 180*phis/pi;

figure(1)
plot(phi_deg, efficiencies)
xlabel('scattering angle, \phi / [deg]')
ylabel('reflectance efficiencies [%]')

x3      = linspace(phi_deg(1), phi_deg(end), 500);
s3      = spline(phi_deg, efficiencies, x3);

figure(2)
plot(x3, s3)
xlabel('scattering angle, \phi / [deg]')
ylabel('reflectance efficiencies [%]')

%
if n == 1
    type = 'gauss';
else
    type = 'super';
end

save(   ['efficiency_' type '.mat'], ...
        'x3',  's3');
%%
s3 = real(s3);
save(   ['efficiency_' type '.mat'], ...
        'x3',  's3');
%%
% figure(3)
% plot(x2 + 90, s2, 'r', x3 + 90, s3, 'b')
% xlabel('scattering angle, \phi / [deg]')
% ylabel('reflectance efficiencies [%]')
% legend('Gaussian', 'super-Gaussian')