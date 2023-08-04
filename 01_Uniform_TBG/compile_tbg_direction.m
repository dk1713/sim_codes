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
dn_g    = 2e-3;

% size
grid_size = 2^11;

dngs_num    = zeros(3, 9);
dngs_ana    = zeros(3, grid_size);

refl_dng_num= zeros(3, 9);
refl_dng_ana= zeros(3, grid_size);

lams_num    = zeros(3, 9);
lams_ana    = zeros(3, grid_size);

dire_lam_num= zeros(3, 9);
dire_lam_ana= zeros(3, grid_size);

thes_num    = zeros(3, 9);
thes_ana    = zeros(3, grid_size);

refl_the_num= zeros(3, 9);
refl_the_ana= zeros(3, grid_size);

phis_deg = [60, 90, 120];
for iter = 1:length(phis_deg)
    phi_deg = phis_deg(iter);
    
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
        '_phi_',    num2str(phi_deg),   ...
        ];
    filename    = ['data/tbg_direction'  model_spec '.mat'];
    load(filename);

    % common parameters
    phi     = phi_deg*pi/180;

    % efficiency studies 
    lambdass = linspace(650, 900, grid_size)*1e-3;

    lams_ana(iter,:) = lambdass*1e3;
    lams_num(iter,:) = lambdas;
    dire_lam_ana(iter,:) = phi_deg*ones(size(lambdass));
    dire_lam_num(iter,:) = phis*180/pi;
%     refl_lam_num(iter,:) = 100*(1 - power_out_lam);
end

%% Figures
figure(4)
plot( ...
    lams_ana(1,:), dire_lam_ana(1,:), 'k', ... 
    lams_ana(2,:), dire_lam_ana(2,:), 'r--', ...
    lams_ana(3,:), dire_lam_ana(3,:), 'b-.', ...
    lams_num(1,:), dire_lam_num(1,:), 'kx', ...
    lams_num(2,:), dire_lam_num(2,:), 'rs', ... 
    lams_num(3,:), dire_lam_num(3,:), 'bo')
xlabel('Input wavelength, \lambda / [nm]')
ylabel('Scattering direction / [deg]')
legend('60^\circ', '90^\circ', '120^\circ', ...
    '60^\circ', '90^\circ', '120^\circ' )

%% write .mat
% save('dataset_3ab.mat', ...
%     'dngs_ana', 'dngs_num', 'refl_dng_ana', 'refl_dng_num', ...
%     'lams_ana', 'lams_num', 'refl_lam_ana', 'refl_lam_num', ...
%     'thes_ana', 'thes_num', 'refl_the_ana', 'refl_the_num')