%% Init
sigma   = 3;
dn      = 5e-3;
n_cl    = 1.4555;
n_co    = 1.4608;
n_eff   = 1.4635;

phi     = 90*pi/180;
theta   = .5*phi;
dn_g    = 1.2e-3;
lam     = 780*1e-3;
Lam     = lam/(n_cl * cos(phi) + n_eff);
K       = 2 * pi./Lam;
w_0     = 2.2719;

%% 1 dn_g efficiency studies
dn_gss  = linspace(0, 5, 100)*1e-3;
kappa   = dn_gss/n_eff;
beta    = 2 * pi * n_eff/lam;

w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;

alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd * kappa.^2 / w_0 / Lam^2; 
alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
alpha_ana   = alpha_ana * exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);

% In dB
% figure(11)
% plot(   dn_gss, 10*alpha_ana*80/log(10), 'b', ... 
%         dn_gs,  -10*log10(power_loss_dng), 'xr')
% xlabel('grating index modulation, {\Delta}n_g')
% ylabel('reflectance at 100 {\mu}m / [dB]')
% legend('analytical', 'numerical' )
% 
% %%
% % In percentage
% figure(12)
% plot(   dn_gss, 100*(1 - 10.^(-alpha_ana*80/log(10))), 'b', ... 
%         dn_gs,  100*(1 - power_loss_dng), 'xr')
% xlabel('grating index modulation, {\Delta}n_g')
% ylabel('reflectance at 100 {\mu}m / [%]')
% legend('analytical', 'numerical' )

%% 2mm
figure(13)
plot(dn_gss, 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10))), 'linewidth', 2);
xlabel('Index modulation, {\Delta}n_g')
ylabel('Reflectance at 3mm grating[%]')
xlim([.2e-3, .8e-3])
yline(4, '--r', 'linewidth', 2);
xline(4.81e-4, '--r', 'linewidth', 2);
set(gca, 'FontSize', 14);

%% fig 3.5 in dB
data_x = dn_gss;
% data_y = 10*alpha_ana*(3e-3)*1e6/log(10);
data_y = 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10)));
save('figure.mat', 'data_x', 'data_y');

%% fig 3.5 in dB
data_x = dn_gss;
% data_y = 10*alpha_ana*(3e-3)*1e6/log(10);
data_y = 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10)));
% save('figure3.5.mat', 'data_x', 'data_y');

%% dB/1mm
figure(22)
plot(   dn_gss, 10*alpha_ana*(1e-3)*1e6/log(10), 'b')
xlabel('\Delta n_g')
ylabel('outcoupling efficiency [dB]')

%% % /10mm
figure(30);
plot(   dn_gss, 100*(1 - 10.^(-alpha_ana*(10e-3)*1e6/log(10))), ...
    'LineWidth', 2, 'markerSize', 15);
xlabel('index modulation, \Delta n_g')
ylabel('diffraction efficiency at 10mm /[%]')

%% % /3mm
figure(40);
plot(   dn_gss, 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10))), ...
    'LineWidth', 2, 'markerSize', 15);
xlabel('index modulation, \Delta n_g')
ylabel('diffraction efficiency at 3mm /[%]')

%% %/10mm for fixed length
figure(31); %dng= 1e-3;
g_dist = linspace(0,20e-3,100);
plot(   g_dist*1e3, 100*(1 - 10.^(-alpha_ana(21)*(g_dist)*1e6/log(10))), ...
    'LineWidth', 2, 'markerSize', 15);
xlabel('grating length/ [mm]')
ylabel('diffraction efficiency at {\Delta}n_g = 1 \times 10^{-3}/[%]')
%% lambdas efficiency studies 
% lambdass = linspace(650, 900, 2001)*1e-3;
% xlabel('grating index modulation, {\Delta}n_g')
% ylabel('reflectance at 10 mm / [%]')

%% 2 lambdas efficiency studies 
% lambdass = linspace(650, 900, 2^10)*1e-3;
% beta    = 2 * pi * n_eff./lambdass;
% kappa   = dn_g/n_eff;
% 
% w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;
% alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam^2; 
% alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
% alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);
% 
% half_angles = 180*asin(.5*lambdass./Lam)/pi;
% 
% [hm_lam, x_hm_lam, fwhm_lam] = find_fwhm(half_angles, 100*(1 - 10.^(-alpha_ana*80/log(10)));
% fprintf('FWHM (theta/2) = %2.4f [deg]\n', fwhm_lam);
% 
% figure(2)
% plot(   half_angles, 100*(1 - 10.^(-alpha_ana*80/log(10))), 'b', ...
%         180*asin(.5*lambdas*1e-3./Lam)/pi, 100*(1 - power_loss_lambda), 'xr')
% xlabel('half angle, \theta /2 / [deg]')
% ylabel('reflectance at 100 {\mu}m / [%]')
% xline(x_hm_lam(1),      '--k',  'HandleVisibility', 'off')
% xline(x_hm_lam(2),      '--k',  'HandleVisibility', 'off')
% yline(hm_lam,          '--k',  'HandleVisibility', 'off')
% legend('analytical', 'numerical')

% %% 2 lambda
% lambdass = linspace(650, 900, 2^12)*1e-3;
% beta    = 2 * pi * n_eff./lambdass;
% kappa   = dn_g/n_eff;
% 
% w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta)^2;
% alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam^2; 
% alpha_ana   = alpha_ana * sin(phi) / (4*cos(theta)^4);
% alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta)^2 - K).^2);
% 
% [hm_lam, x_hm_lam, fwhm_lam] = find_fwhm(lambdass*1e3, 100*(1 - 10.^(-alpha_ana*80/log(10))));
% fprintf('FWHM (lambda) = %2.4f [nm]\n', fwhm_lam);
% 
% figure(2)
% plot(   lambdass*1e3, 100*(1 - 10.^(-alpha_ana*80/log(10))), 'b', ...
%         lambdas, 100*(1 - power_loss_lambda), 'xr')
% xlabel('input wavelength, \lambda /2 / [nm]')
% ylabel('reflectance at 100 {\mu}m / [%]')
% xline(x_hm_lam(1),      '--k',  'HandleVisibility', 'off')
% xline(x_hm_lam(2),      '--k',  'HandleVisibility', 'off')
% yline(hm_lam,          '--k',  'HandleVisibility', 'off')
% legend('analytical', 'numerical')
% 
% %%
% data_x = lambdass*1e3;
% % data_y = 10*alpha_ana*(3e-3)*1e6/log(10);
% data_y = 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10)));
% % save('figure3.6.mat', 'data_x', 'data_y');
% 
% %% 3 tilt angle efficiency studies 
% thetass = linspace(40, 50, 2^10)*pi/180;
% beta    = 2 * pi * n_eff./lam;
% 
% w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*thetass).^2;
% alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam.^2; 
% alpha_ana   = alpha_ana * sin(phi) ./ (4*cos(thetass).^4);
% alpha_ana   = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(thetass).^2 - K).^2);
% 
% [hm_the, x_hm_the, fwhm_the] = find_fwhm(180*thetass/pi, 100*(1 - 10.^(-alpha_ana*80/log(10))));
% fprintf('FWHM (theta) = %2.4f [deg]\n', fwhm_the);
% 
% figure(3)
% plot(   180*thetass/pi, 100*(1 - 10.^(-alpha_ana*80/log(10))), 'b', ...
%         thetas, 100*(1 - power_loss_theta), 'xr')
% xlabel('blazed angle, \theta / [deg]')
% ylabel('reflectance at 100 {\mu}m / [%]')
% xline(x_hm_the(1),      '--k',  'HandleVisibility', 'off')
% xline(x_hm_the(2),      '--k',  'HandleVisibility', 'off')
% yline(hm_the,          '--k',  'HandleVisibility', 'off')
% legend('analytical', 'numerical')
% 
% figure(33)
% plot(   180*thetass/pi, 100*(1 - 10.^(-(cos(2*thetass).^2).*alpha_ana*80/log(10))), 'b')
% xlabel('blazed angle, \theta / [deg]')
% ylabel('reflectance at 100 {\mu}m / [%]')
% 
% %% fig 3.8 - 9 in dB
% the_ana     = 180*thetass/pi;
% % refl_ana_s  = 10*alpha_ana*(3e-3)*1e6/log(10);
% refl_ana_s = 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10)));
% % refl_ana_p  = 10*alpha_ana*(3e-3)*1e6/log(10).*cos(2*thetass).^2;
% refl_ana_p = 100*(1 - 10.^(-alpha_ana*(3e-3)*1e6/log(10).*cos(2*thetass).^2));
% 
% % save('figure3.8.mat', 'the_ana', 'refl_ana_s', 'refl_ana_p')
% 
% %% 4 different scattering angles: 60, 90, 120
% kappa   = dn_g/n_eff;
% 
% % at phi = 60 deg
% phi     = 60*pi/180;
% theta1  = linspace(.5*phi - 10*pi/180, .5*phi + 10*pi/180, 2^12);
% Lam     = lam/(n_cl * cos(phi) + n_eff);
% K       = 2 * pi./Lam;
% 
% w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta1).^2;
% alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam.^2; 
% alpha_ana   = alpha_ana * sin(phi) ./ (4*cos(theta1).^4);
% alpha_ana1  = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta1).^2 - K).^2);
% 
% [hm_the1, x_hm_the1, fwhm_the] = find_fwhm(180*theta1/pi, 100*(1 - 10.^(-alpha_ana1*(1e-3)*1e6/log(10))));
% fprintf('FWHM (phi = 60)  = %2.4f [deg]\n', fwhm_the);
% 
% % at phi = 90 deg
% phi     = 90*pi/180;
% theta2  = linspace(.5*phi - 10*pi/180, .5*phi + 10*pi/180, 2^12);
% Lam     = lam/(n_cl * cos(phi) + n_eff);
% K       = 2 * pi./Lam;
% 
% w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta2).^2;
% alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam.^2; 
% alpha_ana   = alpha_ana * sin(phi) ./ (4*cos(theta2).^4);
% alpha_ana2  = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta2).^2 - K).^2);
% 
% [hm_the2, x_hm_the2, fwhm_the] = find_fwhm(180*theta2/pi, 100*(1 - 10.^(-alpha_ana2*(1e-3)*1e6/log(10))));
% fprintf('FWHM (phi = 90)  = %2.4f [deg]\n', fwhm_the);
% 
% % at phi = 120 deg
% phi     = 120*pi/180;
% theta3  = linspace(.5*phi - 10*pi/180, .5*phi + 10*pi/180, 2^12);
% Lam     = lam/(n_cl * cos(phi) + n_eff);
% K       = 2 * pi./Lam;
% 
% w_th_sqd    = sigma^2 * w_0^2 / (sigma^2 + w_0^2) ./ sin(2*theta3).^2;
% alpha_ana   = pi^2 * sqrt(2*pi) * w_th_sqd .* kappa^2 / w_0 ./ Lam.^2; 
% alpha_ana   = alpha_ana * sin(phi) ./ (4*cos(theta3).^4);
% alpha_ana3  = alpha_ana .* exp(-.5*w_th_sqd .* (2*beta*cos(theta3).^2 - K).^2);
% 
% [hm_the3, x_hm_the3, fwhm_the] = find_fwhm(180*theta3/pi, 100*(1 - 10.^(-alpha_ana3*(1e-3)*1e6/log(10))));
% fprintf('FWHM (phi = 120) = %2.4f [deg]\n', fwhm_the);
% 
% figure(4)
% plot( ...
%     180*theta1/pi, 100*(1 - 10.^(-alpha_ana1*(1e-3)*1e6/log(10))), '--b',   ...
%     180*theta2/pi, 100*(1 - 10.^(-alpha_ana2*(1e-3)*1e6/log(10))), 'k',     ...
%     180*theta3/pi, 100*(1 - 10.^(-alpha_ana3*(1e-3)*1e6/log(10))), '-.r');
% legend('\phi = 60^\circ', '\phi = 90^\circ', '\phi = 120^\circ', 'location', 'SouthEast')
% xlabel('tilt angle, \theta / [deg]')
% ylabel('reflectance at 1 mm / [%]')
% 
% data_x_60   = 180*theta1/pi;
% data_x_90   = 180*theta2/pi;
% data_x_120  = 180*theta3/pi;
% % data_y_60   = 10*alpha_ana1*(3e-3)*1e6/log(10);
% % data_y_90   = 10*alpha_ana2*(3e-3)*1e6/log(10);
% % data_y_120  = 10*alpha_ana3*(3e-3)*1e6/log(10);
% data_y_60   = 100*(1 - 10.^(-alpha_ana1*(3e-3)*1e6/log(10)));
% data_y_90   = 100*(1 - 10.^(-alpha_ana2*(3e-3)*1e6/log(10)));
% data_y_120  = 100*(1 - 10.^(-alpha_ana3*(3e-3)*1e6/log(10)));
% save('figure3b.mat', 'data_x_60', 'data_x_90', 'data_x_120', 'data_y_60', 'data_y_90', 'data_y_120');