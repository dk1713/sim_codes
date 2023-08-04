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
dn_g    = 1e-3;

% size
grid_size = 2^11;

dngs_num    = zeros(3, 9);
dngs_ana    = zeros(3, grid_size);

refl_dng_num= zeros(3, 9);
refl_dng_ana= zeros(3, grid_size);

lams_num    = zeros(3, 9);
lams_ana    = zeros(3, grid_size);

refl_lam_num= zeros(3, 9);
refl_lam_ana= zeros(3, grid_size);

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
    filename    = ['data/tbg_efficiency'  model_spec '.mat'];
    load(filename);

    % common parameters
    phi     = phi_deg*pi/180;
    Lam     = lam/(n_cl * cos(phi) + n_eff);
    K       = 2 * pi./Lam;
    theta   = .5*acos(lam/n_eff/Lam - 1);

    % 1 dn_g efficiency studies
    dn_gss  = linspace(0, 5, grid_size)*1e-3;
    kappa   = dn_gss/n_eff;
    beta    = 2 * pi * n_eff/lam;

    w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;
    alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * kappa.^2 / w_0 / Lam^2; 
    alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
    alpha_ana   = alpha_ana * exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);

    dngs_ana(iter,:) = dn_gss;
    dngs_num(iter,:) = dn_gs;
    refl_dng_ana(iter,:) = 100*(1 - 10.^(-alpha_ana*80/log(10)));
    refl_dng_num(iter,:) = 100*(1 - power_out_dng./power_in);

    % 2 lambdas efficiency studies 
    lambdass = linspace(650, 900, grid_size)*1e-3;
    beta    = 2 * pi * n_eff./lambdass;
    kappa   = dn_g/n_eff;

    w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;
    alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam^2; 
    alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
    alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);

    lams_ana(iter,:) = lambdass*1e3;
    lams_num(iter,:) = lambdas;
    refl_lam_ana(iter,:) = 100*(1 - 10.^(-alpha_ana*80/log(10)));
%     refl_lam_num(iter,:) = 100*(1 - power_out_lam./power_in);
    refl_lam_num(iter,:) = 100*(1 - power_out_lam);
    
    [hm_lam, x_hm_lam, fwhm_lam] = find_fwhm(lambdass*1e3, refl_lam_ana(iter,:));
    fprintf('FWHM (lambda)  = %2.2f     [nm]\n', fwhm_lam);

    % 3 tilt angle efficiency studies 
    thetass = linspace(.5*phi_deg - 10, .5*phi_deg + 10, grid_size)*pi/180;
    beta    = 2 * pi * n_eff./lam;

    w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*thetass).^2;
    alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam.^2; 
    alpha_ana   = alpha_ana * sin(phi) ./ (4*cos(thetass).^4);
    alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(thetass).^2 - K).^2);

    thes_ana(iter,:) = thetass*180/pi;
    thes_num(iter,:) = thetas;
    refl_the_ana(iter,:) = 100*(1 - 10.^(-alpha_ana*80/log(10)));
    refl_the_num(iter,:) = 100*(1 - power_out_the./power_in);
    
    [hm_the, x_hm_the, fwhm_the] = find_fwhm(thes_ana(iter,:), refl_the_ana(iter,:));
    fprintf('FWHM (theta)   = %2.4f     [deg]\n', fwhm_the);
end

%% Figures
figure(1)
plot( ...
    dngs_ana(1,:), refl_dng_ana(1,:), 'k', ... 
    dngs_ana(2,:), refl_dng_ana(2,:), 'r--', ... 
    dngs_num(1,:), refl_dng_num(1,:), 'kx', ...
    dngs_num(2,:), refl_dng_num(2,:), 'rs', ... 
    dngs_num(3,:), refl_dng_num(3,:), 'bo')
xlabel('Index modulation, {\Delta}n_g')
ylabel('Reflectance / [%]')
legend('60^\circ, 120^\circ', '90^\circ', ...
    '60^\circ', '90^\circ', '120^\circ', 'location', 'northwest')

figure(2)
plot( ...
    lams_ana(1,:), refl_lam_ana(1,:), 'k', ... 
    lams_ana(2,:), refl_lam_ana(2,:), 'r--', ...
    lams_ana(3,:), refl_lam_ana(3,:), 'b-.', ...
    lams_num(1,:), refl_lam_num(1,:), 'kx', ...
    lams_num(2,:), refl_lam_num(2,:), 'rs', ... 
    lams_num(3,:), refl_lam_num(3,:), 'bo')
xlabel('Input wavelength, \lambda / [nm]')
ylabel('Reflectance / [%]')
legend('60^\circ', '90^\circ', '120^\circ', ...
    '60^\circ', '90^\circ', '120^\circ' )

figure(3)
plot( ...
    thes_ana(1,:), refl_the_ana(1,:), 'k', ... 
    thes_ana(2,:), refl_the_ana(2,:), 'r--', ... 
    thes_ana(3,:), refl_the_ana(3,:), 'b-.', ... 
    thes_num(1,:), refl_the_num(1,:), 'kx', ...
    thes_num(2,:), refl_the_num(2,:), 'rs', ...
    thes_num(3,:), refl_the_num(3,:), 'bo')
xlabel('Tilt angle, \theta / [deg]')
ylabel('Reflectance / [%]')
legend('60^\circ', '90^\circ', '120^\circ', ...
    '60^\circ', '90^\circ', '120^\circ' )

%% write .mat
% save('dataset_3ab.mat', ...
%     'dngs_ana', 'dngs_num', 'refl_dng_ana', 'refl_dng_num', ...
%     'lams_ana', 'lams_num', 'refl_lam_ana', 'refl_lam_num', ...
%     'thes_ana', 'thes_num', 'refl_the_ana', 'refl_the_num')