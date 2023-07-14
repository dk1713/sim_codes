phis_deg = [60, 90, 120];

lams_num    = zeros(1, 9);
lams_ana    = zeros(1, 2^9);

refl_lam_num= zeros(1, 9);
refl_lam_ana= zeros(1, 2^9);

L       = 120;
H       = 35;
L_g     = 100;
H_g     = 5;
n_cl    = 1.4555;
n_co    = 1.4608;
n_eff   = 1.4640;
sigma   = 3;
dn      = 5e-3;
dn_g    = 1.2e-3;

index   = linspace(n_co, n_co + dn + dn_g, 11);
i_iter  = 11;
% n_1     = index(i_iter);
n_1     = 1.45;

for iter = 1:length(phis_deg)
    phi_deg = phis_deg(iter);
    
    % Load the grating data
%     model_spec = [...
%         '_L_',      num2str(L), ...
%         '_H_',      num2str(H), ...
%         '_L_g_',    num2str(L_g), ...
%         '_H_g_',    num2str(H_g), ...
%         '_n_cl_',   num2str(n_cl), ...
%         '_n_co_',   num2str(n_co), ...
%         '_sigma_',  num2str(sigma), ...
%         '_dn_',     num2str(dn), ...
%         '_phi_',    num2str(phi_deg),   ...
%         '_dng_',    num2str(dn_g), ...
%         ];
%     filename    = ['data/tbg_lambda'  model_spec '.mat'];
%     load(filename);
%     
    lambdass = linspace(650, 900, 2^9)*1e-3;
%     
%     % w_0 extrapolations
    w_0     = interp1(lambdas*1e-3, w_0s, lambdass, 'spline', 'extrap');
%     
%     % phi extrapolations
    phi     = interp1(lambdas*1e-3, phis, lambdass, 'spline', 'extrap');
    
%     figure(10 + iter)
%     plot(lambdass*1e3, sin(phi), lambdas, sin(phis), 'x')
%     yline(sin(phi_deg*pi/180), '--r');
%     xlabel('input wavelength, \lambda / [nm]')
%     ylabel('effective index')
%     ylim([0,1.1]);

    % common parameters
    lam     = 780*1e-3;
    n_eff   = 1.4640;
    Lam     = lam./(n_cl .* cos(phi_deg*pi/180) + n_eff);
    K       = 2 * pi./Lam;
    theta   = .5*acos(lam./n_eff/Lam - 1);
    
    % n_eff extrapolations
%     n_eff   = interp1(lambdas*1e-3, n_effs, lambdass, 'spline', 'extrap');
    
    % Calculating new phi according to propagation constant
    sin_phi = sqrt(1 - ((n_eff - lam./Lam)/n_1).^2);

    % lambdas efficiency studies 
    beta    = 2 * pi * n_1./lambdass;
%     beta    = 2 * pi * n_eff./lambdass;
    kappa   = dn_g./n_eff;

    w_th_sqd    = sigma^2 .* w_0.^2 ./ (sigma^2 + w_0.^2) ./ sin(2*theta)^2;
    alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa.^2 ./ w_0 ./ Lam.^2; 
    alpha_ana   = alpha_ana .* sin(phi) / (4*cos(theta)^4);
%     alpha_ana   = alpha_ana .* sin_phi / (4*cos(theta)^4);
    alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);

    lams_ana(iter,:) = lambdass*1e3;
%     lams_num(iter,:) = lambdas;
    refl_lam_ana(iter,:) = 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10)));
%     refl_lam_ana(iter,:) = 10*alpha_ana*(3e-3)*1e6/log(10);
%     refl_lam_num(iter,:) = 100*(1 - power_loss_lambda);

%     figure(20 + iter)
%     plot(lambdass*1e3, w_0, lambdas, w_0s, 'x')
%     xlabel('input wavelength, \lambda / [nm]')
%     ylabel('guide mode width, w0 / [{\mu}m]')
%     
%     figure(30 + iter)
%     plot(lambdass*1e3, n_eff, lambdas, n_effs, 'x')
%     xlabel('input wavelength, \lambda / [nm]')
%     ylabel('effective index')

end

%% Figures
figure(1)
plot(lambdass*1e3, n_eff, lambdas, n_effs, 'x')
xlabel('input wavelength, \lambda / [nm]')
ylabel('effective index')

figure(10 + i_iter)
plot( ...
    lams_ana(1,:), refl_lam_ana(1,:), 'k', ... 
    lams_ana(2,:), refl_lam_ana(2,:), 'r--', ...
    lams_ana(3,:), refl_lam_ana(3,:), 'b-.')
xlabel('input wavelength, \lambda / [nm]')
ylabel('reflectance / [%]')
legend('60^\circ', '90^\circ', '120^\circ', ...
    '60^\circ', '90^\circ', '120^\circ' )

data_x_60   = lams_ana(1,:);
data_x_90   = lams_ana(2,:);
data_x_120  = lams_ana(3,:);
data_y_60   = refl_lam_ana(1,:);
data_y_90   = refl_lam_ana(2,:);
data_y_120  = refl_lam_ana(3,:);
save('figure3a.mat', 'data_x_60', 'data_x_90', 'data_x_120', 'data_y_60', 'data_y_90', 'data_y_120');