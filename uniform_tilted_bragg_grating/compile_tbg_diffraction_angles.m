% Parameters
L_g     = 100;
L       = L_g + 20;
H       = 35;
H_g     = 5;
n_cl    = 1.4555;
n_co    = 1.4608;
sigma   = 3;
dn      = 5e-3;
dn_g    = 4e-3;

% Load the grating data
model_spec = [...
    '_L_',      num2str(L), ...
    '_H_',      num2str(H), ...
    '_L_g_',    num2str(L_g), ...
    '_H_g_',    num2str(H_g), ...
    '_n_cl_',   num2str(n_cl), ...
    '_n_co_',   num2str(n_co), ...
    '_sigma_',  num2str(sigma), ...
    '_dn_',     num2str(dn), ...
    '_dng_',    num2str(dn_g), ...
    ];
filename    = ['data/tbg_tran_phi'  model_spec '.mat'];
load(filename);

phis_num    = phis_deg*pi/180;
thetas_num  = thetas;

% splining phis
phis_ana    = linspace(phis_num(1), phis_num(end), 2^9);
thetas_ana  = interp1(phis_num, thetas_num, phis_ana, 'spline');

% Common parameters
lam     = 780*1e-3;
Lam     = lam./(n_cl * cos(phis_ana) + n_eff);
K       = 2*pi./Lam;

% Analytic solutions
kappa   = dn_g/n_eff;
beta    = 2 * pi * n_eff/lam;

w_th_sqd    = sigma^2 * w_0^2 ./ (sigma^2 + w_0^2) ./ sin(2*thetas_ana).^2;

alpha_ana   = pi^2 * sqrt(2*pi) .* w_th_sqd .* kappa.^2 / w_0 ./ Lam.^2; 
alpha_ana   = alpha_ana .* sin(phis_ana) ./ (4*cos(thetas_ana).^4);
alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(thetas_ana).^2 - K).^2);

% p-pol
alpha_ana   = alpha_ana .* cos(2*thetas_ana).^2;

refl_DA_ana = 100*(1 - 10.^(-alpha_ana*80/log(10)));
refl_DA_num = 100*(1 - trans);

%% Figures
figure(11)
plot( ...
    phis_ana*180/pi - 90, refl_DA_ana, 'k-',  ... 
    phis_num*180/pi - 90, refl_DA_num, 'bx', 'linewidth', 1, 'markersize', 10)
xlabel('Diffraction angles / [deg]')
ylabel('Reflectance / [%]')
legend('Analitical', 'Numerical', 'location', 'southwest')