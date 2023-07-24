%% Design and test an arbitrary beam generation

phi_deg = 80;               % tilt angle of Gaussian in deg
phi     = phi_deg * pi/180; % tilt angle of Gaussian (90deg = vertical)

d       = 100e-6;           % focal distance
lam     = 0.78e-6;           % wavelength

n_co    = 1.447;
k0      = 2*pi*n_co/lam;
w_0     = 15e-6;            % beam waist
n       = 8;               % order of super Gaussian
outP    = 0.05;

x_tar   = -d /tan(phi);
% defining the range for the grating where .mph file as grating length of
% 100 um, so  x_max - x_min need to be larger.
x_min   = -100e-6;
x_max   = 100e-6;
x  = linspace(x_min, x_max, 5000)';
x_coor = x;

glass_h = 10e-6;

%% Arbitrary field on target location
super_gaussian = @(x, phi, x_tar, w_0, n) ...
    exp(1i*k0*cos(phi) * x) .* exp(-( (x-x_tar)/w_0 ).^( 2*n ));

E_tar   = super_gaussian(x, phi, x_tar, w_0, n);
Ek_tar  = fftshift( fft( fftshift(E_tar) ) );
%%
figure(1)
plot(x, abs(E_tar))
xlabel('x')
ylim([0 1.1])
ylabel('|E|')
title('Field on target')

%% E on the grating
k       = linspace(-pi, pi, length(x) )'./(x(1) - x(2));
ky      = real( sqrt(k0^2 - k.^2) );

Ek_grat = exp(-1i * ky * -(d-glass_h)).*(Ek_tar);
Ez_grat = fftshift(  ifft( fftshift(Ek_grat) ) );

figure(2)
plot(x, abs(Ez_grat))
xline(-60e-6)
xline(60e-6)
xlabel('x')
ylabel('|E|')
title('Field on grating')

%% Parameters for grating
% Remove phase from the tilt
Ez_grat = exp(-1i*k0*cos(phi)*x) .* Ez_grat;

Pz_amp  = abs(Ez_grat).^2;
phase   = unwrap(angle(Ez_grat)); %add 'unwrap' for smoother interpolation!

figure(3)
plot(x, phase)
xlabel('x')
ylabel('phase')
title('desired phase on the grating')

% computing dng_amp depends on pump depletion
F       = griddedInterpolant(x, Pz_amp, 'spline');
fun     = @(x) F(x);
C       = 1/integral(fun, x_min, x_max);

denu    = zeros(size(x));

for i = 1:length(x)
    denu(i) = 1 - outP * C * integral(fun, -100e-6, x(i));
end

dng_amp = sqrt(outP * C * fun(x) ./ denu);

save('phase.mat',       'x',  'phase');
save('dng_amp.mat',     'x',  'dng_amp');

%% Run Comsol
% model = mphload('bragg_grating_design_5.0 (arbitrary field).mph');
% model.param('par2').set('lambda0', [num2str(lam*1e6) '[um]']);
% model.param('par3').set('phi_deg', num2str(phi_deg));
% model.param('default').set('L_g', '120 [um]');
% model.study('std1').run
% 
% % Dimensions
% H       = mphglobal(model, 'H');
% L       = mphglobal(model, 'L');
% L_g     = mphglobal(model, 'L_g');
% phi     = mphglobal(model, 'phi_deg');
% 
% lam     = mphglobal(model, 'lambda0');
% n_clad  = mphglobal(model, 'n_co');
% dn      = mphglobal(model, 'dn');
% 
% x       = linspace(-.5*L, .5*L, 5000);
% xx      = [x; .5*(H - 1)*ones(size(x))];
% Ez      = mphinterp(model, 'ewfd.Ez',     'coord', xx);
% 
% xx      = [x; zeros(size(x))];
% dng_amp = mphinterp(model, 'dng_amp',     'coord', xx);
% max_dng = max(dng_amp);
% 
% y       = linspace(-.5*H, .5*H, 5000);
% yy      = [-.5*L*ones(size(y)); y];
% Px      = mphinterp(model, 'ewfd.Poavx',  'coord', yy);
% p0      = trapz(y, Px);
% 
% yy      = [.5*L*ones(size(y)); y];
% Px      = mphinterp(model, 'ewfd.Poavx',  'coord', yy);
% p1      = trapz(y, Px);
% 
% pow_out = 100*(p0 - p1)/p0;
% 
% %%
% figure(4);
% plot(x, abs(Ez));
% xline(-L_g/2)
% xline(L_g/2)
% xlabel('x')
% ylabel('|E|')
% title('Field on grating')
% 
% save([  'data/bragg_grating',     ...
%         '_dimensions_',         num2str(H), 'x', num2str(L),    ...
%         '_grat_length_',        num2str(L_g),                   ...
%         '_lambda_',             num2str(lam),                   ...
%         '_n_clad_',             num2str(n_clad),                ...
%         '_dn_',                 num2str(dn),                    ...
%         '_scattered_tilt_',     num2str(phi),                   ...
%         '_target_dist_',        num2str(d*1e6),                 ...
%         '_target_waist_',       num2str(w_0*1e6),               ...
%         '_gausssian_order_',    num2str(n),                     ...
%         '_power_ratio_',        num2str(outP),                  ...
%         '.mat'], 'x', 'Ez', 'pow_out', 'max_dng');