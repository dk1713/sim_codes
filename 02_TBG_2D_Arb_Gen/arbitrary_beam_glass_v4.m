%% Design and test an arbitrary beam generation
% Indices
n_air   = 1;
n_clad  = 1.537;
n_core  = 1.548;
% n_under = 1.547;
dn_g    = 3e-3;

% Heights of the layers [m]
h_clad  = 15e-6;
h_core  = 5e-6;

% Parameters of the Gaussians
% Positions of the Gaussian addressing [m]
d       = 60 * 1e-6;
r_1     = [-30*1e-6, d];
r_2     = [30*1e-6, d];
% Common features of the Gaussians
w       = 2e-6;
n       = 1;
lambda  = 780e-9;
outP    = 0.2;
% Features on the grating
w_g     = 10e-6;

% wavenumbers
k0      = 2*pi*n_air/lambda;
k0_clad = 2*pi*n_clad/lambda;
k0_core = 2*pi*n_core/lambda;

%% Define beam in target plane
% tilt angle of Gaussians (0 rad = vertical)
phis    = atan([r_1(1)/(r_1(2) + h_clad + h_core/2),    ...
                r_2(1)/(r_2(2) + h_clad + h_core/2)]);

% Calculate the new phi's
r_xs    = [r_1(1), r_2(1)];
r_dist  = [0, 0];
while and(mean(abs((r_xs - r_dist)./r_xs)) > 1e-6, abs(r_xs) > abs(r_dist))
    phis    = phis*(1 + 1e-6);
    phi2    = asin(n_air*sin(phis)/n_clad);
    phi3    = asin(n_clad*sin(phi2)/n_core);
    r_dist  = .5*h_core*tan(phi3) + h_clad*tan(phi2) + d*tan(phis);
end

% defining the range for the grating where .mph file as grating length of
% 100 um, so  x_max - x_min need to be larger.
x_min   = -100e-6;
x_max   = 100e-6;
x       = linspace(x_min, x_max, 2^11)';
x_ori   = x;

%% Arbitrary field on target location
gaussian_air = @(x, phi, x_tar, w_0, n) ...
    exp(-1i * k0 * (  sin(phi)*x  )) .* exp(-( (x-x_tar)/w_0 ).^( 2*n ));

num_add     = 2;
x_tar       = zeros(1, num_add);
x_tar(1)    = r_1(1);
x_tar(2)    = r_2(1);

Ez_ori  = 0;
for ii = 1:num_add
    Ez_ori = Ez_ori + gaussian_air(x, phis(ii), x_tar(ii), w, n);
end

Ek_tar  = fftshift( fft( fftshift(Ez_ori) ) );
%% Checking the multi-addressing Gaussians
figure(1)
plot(x, abs(Ez_ori))
xlabel('x')
ylim([0 1.1])
ylabel('|E|')
title('Field on target')

%% Propagation onto the grating plane
k       = linspace(-pi, pi, length(x) )'./(x(1) - x(2));
ky_air  = real( sqrt(k0^2 - k.^2) );
ky_clad = real( sqrt(k0_clad^2 - k.^2) );
ky_core = real( sqrt(k0_core^2 - k.^2) );

phase_air   = asin(k/(n_air*k0));
phase_clad  = asin(k/(n_clad*k0));
phase_core  = asin(k/(n_core*k0));

Ek_grat = exp(-1i*ky_air*-d).*Ek_tar;
Ek_grat = apply_Fresnel(Ek_grat, phase_air, phase_clad, n_air, n_clad, 's');
Ek_grat = exp(-1i*ky_clad*-h_clad).*Ek_grat;
Ek_grat = apply_Fresnel(Ek_grat, phase_clad, phase_core, n_clad, n_core, 's');
Ek_grat = exp(-1i*ky_core*-h_core/2).*Ek_grat;

Ez_grat  = fftshift(  ifft( fftshift(Ek_grat) )  );

figure(2)
plot(x, abs(Ez_grat))
xline(-50e-6)
xline(50e-6)
xlabel('x')
ylabel('|E|')
title('Field on grating')

%% Parameters for grating
% Remove phase from the tilt
% psi     = asin(n_clad*asin(phi/n_clad)/n_core);
Ez_grat = exp(1i*k0_clad*sin(mean(phis))*x) .* Ez_grat;

Pz_amp  = abs(Ez_grat).^2;
phase   = unwrap(angle(Ez_grat)); %add 'unwrap' for smoother interpolation!

% computing dng_amp depends on pump depletion
F       = griddedInterpolant(x, Pz_amp, 'spline');
fun     = @(x) F(x);
C       = 1/integral(fun, x_min, x_max);

denu    = zeros(size(x));

for i = 1:length(x)
    denu(i) = 1 - outP * C * integral(fun, -100e-6, x(i));
end

dng_amp = sqrt(outP * C * fun(x) ./ denu);

%%
figure(3)
plot(x, phase)
xlabel('x')
ylabel('phase')
title('desired phase on the grating')
xline(-50e-6)
xline(50e-6)

figure(4)
plot(x, dng_amp)
xlabel('x')
ylabel('grating strength')
title('desired grating strength')
xline(-50e-6)
xline(50e-6)

save('phase.mat',       'x',  'phase');
save('dng_amp.mat',     'x',  'dng_amp');

%% Run Comsol
model = mphload('grating_design_arbitrary_v5.3.mph');
model.param('par2').set('phi', [num2str(90 + mean(phis)*180/pi), ' [deg]']);
model.study('std1').run

% Dimensions
H       = mphglobal(model, 'H');
L       = mphglobal(model, 'L');
L_g     = mphglobal(model, 'L_g');

lam     = mphglobal(model, 'lambda0');
n_clad  = mphglobal(model, 'n_cl');
dn      = mphglobal(model, 'dn');

x       = linspace(-.5*L, .5*L, 5000);
xx      = [x; (.5*H - 3)*ones(size(x))];
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
figure(5);
plot(x, abs(Ez));
xline(-L_g/2)
xline(L_g/2)
xlabel('x')
ylabel('|E|')
title('Field on measurement line')

% save([  'data/grating_design',     ...
%         '_dimensions_',         num2str(H), 'x', num2str(L),    ...
%         '_grat_length_',        num2str(L_g),                   ...
%         '_lambda_',             num2str(lam),                   ...
%         '_n_clad_',             num2str(n_clad),                ...
%         '_dn_',                 num2str(dn),                    ...
%         '_target_dist_',        num2str(d*1e6),                     ...
%         '_target_x1_',          num2str(r_1(1)*1e6),            ...
%         '_target_x2_',          num2str(r_2(1)*1e6),            ...
%         '_target_waist_',       num2str(w*1e6),                 ...
%         '_gausssian_order_',    num2str(n),                     ...
%         '_power_ratio_',        num2str(outP),                  ...
%         '.mat'],    'x', 'x_ori', 'Ez', 'Ez_ori', 'pow_out', 'max_dng', ...
%                     'x0', 'y0', 'Ez_wg');