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
phi_degs = [60, 90, 120];

sigs_ana        = zeros(length(phi_degs), grid_size);
sigs_num        = zeros(length(phi_degs), 10);
refl_sig_ana    = sigs_ana;
refl_sig_num    = sigs_num;

for ite = 1:length(phi_degs)
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
        '_phi_',    num2str(phi_degs(ite)),   ...
        ];
    filename    = ['data/tbg_sigma'  model_spec '.mat'];
    load(filename);

    % common parameters
    phi     = phi_degs(ite)*pi/180;
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

    sigs_ana(ite, :)        = sigmass;
    sigs_num(ite, :)        = sigmas;
    refl_sig_ana(ite, :)    = 100*(1 - 10.^(-alpha_ana*80/log(10)));
    refl_sig_num(ite, :)    = 100*(1 - power_out_sig);
end

%% Figures
figure(1); clf;
plot(sigmas, n_effs, 'x-', 'LineWidth', 2);
xlabel('sigma, {\sigma} / [um]');
ylabel('effective index, n_{eff}');
set(gca, 'FontSize', 16);

figure(2); clf;
plot(sigmas, w_0s, 'x-', 'LineWidth', 2);
xlabel('sigma, {\sigma} / [um]');
ylabel('mode width, w_0 / [um]');
set(gca, 'FontSize', 16);

figure(3); clf;
plot(sigmas, errors, 'x-', 'LineWidth', 2);
xlabel('sigma, {\sigma} / [um]');
ylabel('relative squred error');
set(gca, 'FontSize', 16);

figure(4); clf;
plot(...
    sigs_ana(1, :), refl_sig_ana(1, :), '-', ...
    sigs_ana(2, :), refl_sig_ana(2, :), '-', ...
    sigs_ana(3, :), refl_sig_ana(3, :), '--', ...
    sigs_num(1, :), refl_sig_num(1, :), 'x', ...
    sigs_num(2, :), refl_sig_num(2, :), '+', ...
    sigs_num(3, :), refl_sig_num(3, :), 'o', 'LineWidth', 2, 'markerSize', 15);
xlabel('sigma, {\sigma} / [um]');
ylabel('Diffraction efficiency / [%]');
colororder(["#8040E6";"#1AA640";"#E68000"])
legend(...
    '60^\circ', '90^\circ', '120^\circ', ...
    '60^\circ', '90^\circ', '120^\circ', 'Location', 'northwest');

%% Figure for the paper
% figure 5: o(sig^2 w0 / (w0^2 + sig^2) )
% figure 6: o(1/n_eff^2)

ord_w0 = sigmas.^2.*w_0s./(sigmas.^2 + w_0s.^2);
figure(5); clf;
plot(sigmas, ord_w0, 'x-', 'LineWidth', 2);
xlabel('sigma, {\sigma} / [um]');
ylabel('o(\sigma^2 w_0/(\sigma^2 + w_0^2) / [1/um]');
set(gca, 'FontSize', 16);

ord_n_eff = 1./n_effs.^2;
figure(6); clf;
plot(sigmas, ord_n_eff, 'x-', 'LineWidth', 2);
xlabel('sigma, {\sigma} / [um]');
ylabel('o(1/n_{eff}^2)');
set(gca, 'FontSize', 16);

%% write .mat
% save('dataset_3ab.mat', ...
%     'dngs_ana', 'dngs_num', 'refl_dng_ana', 'refl_dng_num', ...
%     'lams_ana', 'lams_num', 'refl_lam_ana', 'refl_lam_num', ...
%     'thes_ana', 'thes_num', 'refl_the_ana', 'refl_the_num')