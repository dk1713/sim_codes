%% Init
L       = 120;
H       = 35;
L_g     = 100;
H_g     = 5;
n_cl    = 1.4555;
n_co    = 1.4608;
sigma   = 3;
dn      = 5e-3;
lam     = 780*1e-3;

% Changing for dataset depends on dng
dn_g    = 3e-3;

% size
grid_size = 2^11;
phi_deg = 90;
    
%% Load the grating data
model_spec = [...
    '_L_',      num2str(L), ...
    '_H_',      num2str(H), ...
    '_L_g_',    num2str(L_g), ...
    '_H_g_',    num2str(H_g), ...
    '_n_cl_',   num2str(n_cl), ...
    '_n_co_',   num2str(n_co), ...
    '_sigma_',  num2str(sigma), ...
    '_dn_',     num2str(dn), ...
    '_phi_',    num2str(phi_deg),   ...
    ];
filename    = ['data/tbg_sigma'  model_spec '.mat'];
load(filename);

% common parameters
phi     = phi_deg*pi/180;
Lam     = lam/(n_cl * cos(phi) + n_eff);
K       = 2 * pi./Lam;
theta   = .5*acos(lam/n_eff/Lam - 1);

% 1 dn_g efficiency studies
sigmass = linspace(0, 6, grid_size);
kappa   = dn_g/n_eff;
beta    = 2 * pi * n_eff/lam;

w_th_sqd    = sigmass.^2 * w_0^2 ./ (sigmass.^2 + w_0^2) ./ sin(2*theta)^2;
alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * kappa.^2 / w_0 / Lam^2; 
alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);

sigs_ana        = sigmass;
sigs_num        = sigmas;
refl_sig_ana    = 100*(1 - 10.^(-alpha_ana*80/log(10)));
refl_sig_num    = 100*(1 - power_out_sig);

%% Figures
figure(1); clf;
plot(sigmas, n_effs, 'x-');
yline(n_eff);
xlabel('sigma, {\sigma} / [um]');
ylabel('effective index, n_{eff}');

figure(2); clf;
plot(sigmas, w_0s, 'x-');
yline(w_0);
xlabel('sigma, {\sigma} / [um]');
ylabel('waveguide width, w_0 / [um]');

figure(3); clf;
plot(sigmas, errors, 'x-')
xlabel('sigma, {\sigma} / [um]');
ylabel('relative squred error');

figure(4); clf;
plot(...
    sigs_ana, refl_sig_ana, '-', ...
    sigs_num, refl_sig_num, 'x')
xlabel('sigma, {\sigma} / [um]');
ylabel('Reflectance / [%]');

%% write .mat
% save('dataset_3ab.mat', ...
%     'dngs_ana', 'dngs_num', 'refl_dng_ana', 'refl_dng_num', ...
%     'lams_ana', 'lams_num', 'refl_lam_ana', 'refl_lam_num', ...
%     'thes_ana', 'thes_num', 'refl_the_ana', 'refl_the_num')