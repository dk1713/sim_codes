%% Load the grating data
L       = 120;
H       = 35;
L_g     = 100;
H_g     = 5;
n_cl    = 1.4555;
n_co    = 1.4608;
sigma   = 3;
dn      = 5e-3;

model_spec = [...
    '_L_',      num2str(L), ...
    '_H_',      num2str(H), ...
    '_L_g_',    num2str(L_g), ...
    '_H_g_',    num2str(H_g), ...
    '_n_cl_',   num2str(n_cl), ...
    '_n_co_',   num2str(n_co), ...
    '_sigma_',  num2str(sigma), ...
    '_dn_',     num2str(dn), ...
    ];
filename    = ['data/tbg_powerflow'  model_spec '.mat'];
load(filename);

dn_g    = 4e-3;
lam     = 780*1e-3;

%% 60deg
phi     = 60*pi/180;
theta   = 30*pi/180;
Lam     = lam/(n_cl * cos(phi) + n_eff);
K       = 2 * pi./Lam;

dn_gss  = linspace(0, 5, 100)*1e-3;
kappa   = dn_gss/n_eff;
beta    = 2 * pi * n_eff/lam;

w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;

alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * kappa.^2 / w_0 / Lam^2; 
alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
alpha_ana   = alpha_ana * exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);
alpha_60    = 100*(1 - 10.^(-alpha_ana*(10e-3)*1e6/log(10)));

%% 90deg
phi     = 90*pi/180;
theta   = 45*pi/180;
Lam     = lam/(n_cl * cos(phi) + n_eff);
K       = 2 * pi./Lam;

dn_gss  = linspace(0, 5, 100)*1e-3;
kappa   = dn_gss/n_eff;
beta    = 2 * pi * n_eff/lam;

w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;

alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * kappa.^2 / w_0 / Lam^2; 
alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
alpha_ana   = alpha_ana * exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);
alpha_90    = 100*(1 - 10.^(-alpha_ana*(10e-3)*1e6/log(10)));

%% 120deg
phi     = 120*pi/180;
theta   = 60*pi/180;
Lam     = lam/(n_cl * cos(phi) + n_eff);
K       = 2 * pi./Lam;

dn_gss  = linspace(0, 5, 100)*1e-3;
kappa   = dn_gss/n_eff;
beta    = 2 * pi * n_eff/lam;

w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;

alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * kappa.^2 / w_0 / Lam^2; 
alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
alpha_ana   = alpha_ana * exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);
alpha_120   = 100*(1 - 10.^(-alpha_ana*(10e-3)*1e6/log(10)));

% 10m
figure(13)
plot(   dn_gss, alpha_60,   'k', ...
        dn_gss, alpha_90,   'r--')
xlabel('Index modulation, {\Delta}n_g')
ylabel('Reflectance [%]')
legend('60^\circ, 120^\circ', '90^\circ')