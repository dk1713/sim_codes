%% Fixed lengths calculations for uniform grating.
% Looking into the optimum efficiency for varying scattering direction,
% phi. This script is created to look at uniform grating dataset around 2
% parameter lengths:
%       1. fixed grating length, and
%       2. fixed width of the beam.
% and while maintaining the optimum grating period and tilt angles
% (constant). The dimension for the length is [um] only because COMSOL
% simulations before was so.

%% Init
lam     = 780e-3;
sigma   = 3;
n_cl    = 1.4555;
n_co    = 1.4608;
dn      = 5e-3;
dn_g    = 1.2e-3;
w_0     = 2.2719;

n_eff   = 1.4640; % This is from COMSOL simulation

%% Computing the efficiency
len_g   = (2e-3)*1e6; % 2mm
phis_deg= linspace(50, 80, 2^9);
phis    = phis_deg*pi/180;

theta   = .5*phis; 
% theta is approx estimate but should be close enough according to the
% theory. (Yoshino, Posner, Ko)
Lam     = lam./(n_cl * cos(phis) + n_eff);
K       = 2 * pi./Lam;

kappa       = dn_g/n_eff;
beta        = 2 * pi * n_eff/lam;
w_the_sqd   = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta).^2;

alpha_ana = pi^2 * sqrt(2*pi) * w_the_sqd * kappa.^2 / w_0 / Lam.^2; 
alpha_ana = alpha_ana .* sin(phis) ./ (4*cos(theta).^4);
alpha_ana = alpha_ana .* exp(-.5*w_the_sqd .* (2*beta*cos(theta).^2 - K).^2);

%% 1. Fixed grating length
eff_fixed_lg_db  = 10*alpha_ana*len_g/log(10);
eff_fixed_lg_per = 100*(1 - 10.^(-alpha_ana*len_g/log(10)));

% In dB
figure(11); clf;
plot( phis_deg, eff_fixed_lg_db);
xlabel('grating index modulation, {\Delta}n_g');
ylabel('reflectance at 100 {\mu}m / [dB]');

% In percentage
figure(12); clf;
plot( phis_deg, eff_fixed_lg_per);
xlabel('grating index modulation, {\Delta}n_g');
ylabel('reflectance at 100 {\mu}m / [%]');

%% 2. Fixed beam width
beam_wid = len_g;
len_g = beam_wid./sin(phis);

eff_fixed_w_db   = 10*alpha_ana.*len_g/log(10);
eff_fixed_w_per  = 100*(1 - 10.^(-alpha_ana.*len_g/log(10)));

% In dB
figure(21); clf;
plot( phis_deg, eff_fixed_w_db);
xlabel('scattering direction, \phi / [deg]');
ylabel('reflectance at 100 {\mu}m / [dB]');

% In percentage
figure(22); clf;
plot( phis_deg, eff_fixed_w_per);
xlabel('scattering direction, \phi / [deg]');
ylabel('reflectance at 100 {\mu}m / [%]');

%% save data
% save('in_plane_efficiencies.mat', 'phis', ...
%     'eff_fixed_lg_db', 'eff_fixed_lg_per', ...
%     'eff_fixed_w_db',  'eff_fixed_w_per')