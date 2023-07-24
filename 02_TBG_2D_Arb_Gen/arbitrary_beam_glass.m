%% Design and test an arbitrary beam generation
n_co    = 1.447;
% n_co    = 1.4459;
% n_cl    = 1.4454;
w_0     = 15e-6;            % beam waist
n       = 10;               % order of super Gaussian
outP    = 0.1;

%% Analysis of Fresnel's coefficients (s-pol only)
phi     = linspace(-pi/2, pi/2, 200);

sin_g   = @(phi) sin(phi)./n_co;
cos_g   = @(phi) sqrt(1 - sin_g(phi).^2);
ref_s   = @(phi) (  (n_co*cos_g(phi) - cos(phi)) ... 
                    ./ (n_co*cos_g(phi) + cos(phi)) );

R_s     = ref_s(phi).^2;
T_s     = 1 - R_s;

phi_lim_deg = 30;
phi_lim     = pi*phi_lim_deg/180;
T_value     = 1 - ref_s(phi_lim)^2;

figure(1)
plot(phi*180/pi, R_s, 'b', phi*180/pi, T_s, 'k')
xlabel('Scattered angle [deg]')
ylabel('Fresnel coefficients [1]')
title('Fresnel coefficients (s-pol)')
legend('reflection', 'transmission')
xline(-90,              'k',    'HandleVisibility', 'off')
xline(90,               'k',    'HandleVisibility', 'off')
xline(-phi_lim_deg,     '--k',  'HandleVisibility', 'off')
xline(phi_lim_deg,      '--k',  'HandleVisibility', 'off')
yline(T_value,          '--r',  'HandleVisibility', 'off')

disp('Fresnel Coefficients Analysis')
disp(['  phi bandwidth     [deg]:  ' ...
    '[' num2str(-phi_lim_deg) ', ' num2str(phi_lim_deg) ']'])
disp(['  upper T_s           [%]:  ' num2str(100*T_value)])
disp(' ')

%% Define parameters in air
% choose angle such that it minimises the reflection as much as possible.
phi_deg = -20;              % tilt angle of Gaussian in deg
phi     = phi_deg * pi/180; % tilt angle of Gaussian (0 rad = vertical)

h       = 30e-6;            % height of glass boundary
d       = 80e-6;            % focal distance
lam     = 0.78e-6;          % wavelength
k0_a    = 2*pi/lam;         % k in air
k0_g    = 2*pi*n_co/lam;    % k in grating

y_tar   = d;
x_tar   = d * tan(phi);
disp(['  target coordinates [um]:  ' ...
    '[' num2str(1e6*x_tar) ', ' num2str(1e6*y_tar) ']'])

% defining the range for the grating where .mph file as grating length of
% 100 um, so  x_max - x_min need to be larger.
x_min   = -100e-6;
x_max   = 100e-6;
x       = linspace(x_min, x_max, 5000)';
x_ori   = x;

%% Arbitrary field on target location
super_gaussian = @(x, phi, x_tar, w_0, n) ...
    exp(-1i * k0_a * sin(phi) * x) .* exp(-( (x-x_tar)/w_0 ).^( 2*n ));

Ez_ori  = super_gaussian(x, phi, x_tar, w_0, n);
Ek_tar  = fftshift( fft( fftshift(Ez_ori) ) );
%%
figure(2)
plot(x*1e6, abs(Ez_ori))
xlabel('x [\mu m]')
ylim([0 1.1])
ylabel('|E| [Vm^{-1}]')
title('Field on target')
xline(-60,    '--k', 'HandleVisibility','off')
xline(60,     '--k', 'HandleVisibility','off')

%% normal propagation
k       = linspace(-pi, pi, length(x) )'./(x(1) - x(2));
ky_a    = real( sqrt(k0_a^2 - k.^2) );
ky_g    = real( sqrt(k0_g^2 - k.^2) );

y       = linspace(-y_tar, 0, 400);
Ek_pro  = exp(-1i*ky_a*y).*(Ek_tar*ones(size(y)));

y2      = linspace(-h, 0, 400);
Ek_pro2 = exp(-1i*ky_g*y2).*(Ek_pro(:,1)*ones(size(y2)));

Ek_pro  = [Ek_pro2 Ek_pro];
Ez_pro  = fftshift(  ifft( fftshift(Ek_pro,1) ),1  );

figure(3)
pcolor(x*1e6, [y2+h y+d+h]*1e6, abs(Ez_pro)')
xlabel('x [\mu m]')
ylabel('y [\mu m]')
shading flat
colorbar
axis equal
yline(d + h,    '--k', 'HandleVisibility','off')
yline(h,    	'--r', 'HandleVisibility','off')
yline(0,        '--k', 'HandleVisibility','off')

%% on the grating
Ek_grat = exp(-1i*ky_g*-h -1i*ky_a*-y_tar).*Ek_tar;
Ez_grat = fftshift(  ifft( fftshift(Ek_grat) ) );

figure(4)
plot(x*1e6, abs(Ez_grat))
xline(-60)
xline(60)
xlabel('x [\mu m]')
ylabel('|E| [Vm^{-1}]')
title('Field on grating')

% Parameters for grating
% Remove phase from the tilt
Ez_grat = exp(1i*k0_g*sin(phi)*x) .* Ez_grat;

Pz_amp  = abs(Ez_grat).^2;
phase   = unwrap(angle(Ez_grat)); %add 'unwrap' for smoother interpolation!

figure(5)
plot(x*1e6, phase)
xlabel('x [\mu m]')
ylabel('phase [rad]')
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
model = mphload('bragg_grating_design_5.0 (arbitrary field).mph');
model.param('par3').set('phi_deg', num2str(90 + phi_deg));
model.param('default').set('L_g', '120 [um]');
model.study('std1').run

% Dimensions
H       = mphglobal(model, 'H');
L       = mphglobal(model, 'L');
L_g     = mphglobal(model, 'L_g');
phi     = mphglobal(model, 'phi_deg');

lam     = mphglobal(model, 'lambda0');
n_clad  = mphglobal(model, 'n_co');
dn      = mphglobal(model, 'dn');

x       = linspace(-.5*L, .5*L, 5000);
xx      = [x; .5*(H - 1)*ones(size(x))];
Ez      = mphinterp(model, 'ewfd.Ez',     'coord', xx);

% picture of waveguide
x0       = linspace(-.5*L, .5*L, 2000);
y0       = linspace(-.5*H, .5*H, 2000);
[x1, y1] = meshgrid(x0, y0);
xx       = [x1(:), y1(:)]';
Ez_wg    = mphinterp(model, 'ewfd.Ez',     'coord', xx);
Ez_wg    = reshape(Ez_wg, length(x0), length(y0));

% the efficiency parameters
xx      = [x; zeros(size(x))];
dng_amp = mphinterp(model, 'dng_amp',     'coord', xx);
max_dng = max(dng_amp);

y       = linspace(-.5*H, .5*H, 5000);
yy      = [-.5*L*ones(size(y)); y];
Px      = mphinterp(model, 'ewfd.Poavx',  'coord', yy);
p0      = trapz(y, Px);

yy      = [.5*L*ones(size(y)); y];
Px      = mphinterp(model, 'ewfd.Poavx',  'coord', yy);
p1      = trapz(y, Px);

pow_out = 100*(p0 - p1)/p0;

%%
figure(6);
plot(x, abs(Ez));
xline(-L_g/2)
xline(L_g/2)
xlabel('x [\mu m]')
ylabel('|E| [Vm^{-1}]')
title('Field on measurement line')

save([  'data/bragg_grating',     ...
        '_dimensions_',         num2str(H), 'x', num2str(L),    ...
        '_grat_length_',        num2str(L_g),                   ...
        '_lambda_',             num2str(lam),                   ...
        '_n_clad_',             num2str(n_clad),                ...
        '_dn_',                 num2str(dn),                    ...
        '_scattered_tilt_',     num2str(phi),                   ...
        '_target_dist_',        num2str(d*1e6),                 ...
        '_target_waist_',       num2str(w_0*1e6),               ...
        '_gausssian_order_',    num2str(n),                     ...
        '_power_ratio_',        num2str(outP),                  ...
        '.mat'],    'x', 'x_ori', 'Ez', 'Ez_ori', 'Ez_grat', 'pow_out', 'max_dng', ...
                    'x0', 'y0', 'Ez_wg');