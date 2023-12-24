%% Further calculations for uniform grating.
% Looking into efficiency with grating length

%% Init
lam     = 780e-3;
sigma   = 3;
n_cl    = 1.4555;
n_co    = 1.4608;
dn      = 5e-3;
w_0     = 2.2719;

n_eff   = 1.4640; % This is from COMSOL simulation

%% Computing the efficiency
phi_deg = 90;
phi     = phi_deg*pi/180;

theta   = .5*phi; 
Lam     = lam./(n_cl * cos(phi) + n_eff);
K       = 2 * pi./Lam;

% dng = 1.2e-3
dn_g    = 1.2e-3;
kappa       = dn_g/n_eff;
beta        = 2 * pi * n_eff/lam;
w_the_sqd   = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta).^2;

alpha_ana = pi^2 * sqrt(2*pi) * w_the_sqd * kappa.^2 / w_0 / Lam.^2; 
alpha_ana = alpha_ana .* sin(phi) ./ (4*cos(theta).^4);
alpha_ana1= alpha_ana .* exp(-.5*w_the_sqd .* (2*beta*cos(theta).^2 - K).^2);

% dng = 5.1e-3
dn_g    = 5.1e-3;
kappa       = dn_g/n_eff;
beta        = 2 * pi * n_eff/lam;
w_the_sqd   = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta).^2;

alpha_ana = pi^2 * sqrt(2*pi) * w_the_sqd * kappa.^2 / w_0 / Lam.^2; 
alpha_ana = alpha_ana .* sin(phi) ./ (4*cos(theta).^4);
alpha_ana2= alpha_ana .* exp(-.5*w_the_sqd .* (2*beta*cos(theta).^2 - K).^2);

% dng = 5.0e-4
dn_g    = 5.0e-4;
kappa       = dn_g/n_eff;
beta        = 2 * pi * n_eff/lam;
w_the_sqd   = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta).^2;

alpha_ana = pi^2 * sqrt(2*pi) * w_the_sqd * kappa.^2 / w_0 / Lam.^2; 
alpha_ana = alpha_ana .* sin(phi) ./ (4*cos(theta).^4);
alpha_ana3= alpha_ana .* exp(-.5*w_the_sqd .* (2*beta*cos(theta).^2 - K).^2);

%% figure
lens  = linspace(0, 10, 2^9)*1e3; % [mm]
eff_case1 = 100*(1 - 10.^(-alpha_ana1*lens/log(10)));
eff_case2 = 100*(1 - 10.^(-alpha_ana2*lens/log(10)));
eff_case3 = 100*(1 - 10.^(-alpha_ana3*lens/log(10)));

lens = lens*1e-3;
figure(12); clf;
plot( ...
    lens, eff_case1, ...
    lens, eff_case2, ...
    lens, eff_case3, 'LineWidth', 2 ...
    );
xlabel('grating length / [mm]');
ylabel('reflectance / [%]');
legend( ...
    '{\Delta}n_{ac} = 1.2 \times 10^{-3}', ...
    '{\Delta}n_{ac} = 5.1 \times 10^{-3}', ...
    '{\Delta}n_{ac} = 5.0 \times 10^{-4}');
set(gca, 'FontSize', 14);

%% save data
save('length_efficiencies.mat', 'lens', ...
    'eff_case1', 'eff_case2', 'eff_case3')