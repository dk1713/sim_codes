%% free space propagation and Gaussian beam fitting
% load electric field above waveguide from Comsol simulation
% Indices [1]
n_air   = 1;
n_clad  = 1.4555;
n_core  = 1.4608;
n_eff   = 1.4635;
dn      = 5e-3;

% Heights of the layers [m]
h_clad  = 15e-6;
h_core  = 5e-6;

% Target distance and beam waist
lam     = 780e-9;
r       = [-3, 3];
d       = 100;
w_0     = 2;
n       = 1;
outP    = 0.2;

h       = (h_clad + .5*h_core)*1e6;

% Dimensions
L       = 120;
H       = 30;
L_g     = 100;

load([  'data/grating_design',     ...
        '_dimensions_',         num2str(H), 'x', num2str(L),    ...
        '_grat_length_',        num2str(L_g),                   ...
        '_lambda_',             num2str(lam * 1e6),             ...
        '_n_clad_',             num2str(n_clad),                ...
        '_dn_',                 num2str(dn),                    ...
        '_target_dist_',        num2str(d),                     ...
        '_target_x1_',          num2str(r(1)),              ...
        '_target_x2_',          num2str(r(2)),              ...
        '_target_waist_',       num2str(w_0),                   ...
        '_gausssian_order_',    num2str(n),                     ...
        '_power_ratio_',        num2str(outP),                  ...
        '.mat']);

%% Initialisation
x       = x * 1e-6;
k0      = 2*pi * n_air / lam;
k0_clad = 2*pi * n_clad / lam;
k0_core = 2*pi * n_core / lam;

y_shift = (.5*H - 3)*1e-6;

% Defining the element size
N       = length(x);
dx      = x(2) - x(1);

% Defining the Fourice space
k       = linspace(-pi,pi,N)'/dx;
Ek      = fftshift( fft( fftshift(Ez) ) );

%% Filter the spectrum to get the side-scattered light only
lower_bound = -1*k0;
upper_bound = 1*k0;

index   = find(k>lower_bound & k<upper_bound);
Ek_fil  = zeros(size(Ek));
Ek_fil(index)   = Ek(index);

E_fil   = fftshift( ifft( fftshift(Ek_fil) ) );

%% Padding the x-axis
% now pad the field with zeros in real space (double no. of grid points)

N_fac   = 3;                % factor by which to increase the grid
N2      = N_fac * N;        % new size

k2      =   linspace(-pi,pi,N2)'/dx;
x2      =   ( (-.5*N2:.5*N2-1) + (.5*N+1) )*dx + x(1);

Ez_fil2  = zeros(N2,1);
Ez_fil2(.5*(N2-N) + (0:N-1)) = E_fil;
Ek_fil2 = fftshift( fft( fftshift(Ez_fil2) ) );
%% Electric field in real and fourier spaces
figure(11)
plot(   x,  abs(Ez),        'b', ...
        x2, abs(Ez_fil2),   '--r' )
xlabel('x')
ylabel('|E|')
legend('original', 'filtered')

figure(12)
plot(   k,  abs(Ek),        'b', ...
        k2, abs(Ek_fil2),   '--r')
xlabel('k')
ylabel('|Ek|')
legend('original', 'filtered')
hold on
xline(lower_bound,  '--k','HandleVisibility','off')
xline(upper_bound,  '--k','HandleVisibility','off')
hold off

%% Free Space propagation in the y-axis direction
ky      = real( sqrt(k0^2 - k2.^2) );
ky_clad = real( sqrt(k0_clad^2 - k2.^2) );

phase_air   = asin(k2/(n_air*k0));
phase_clad  = asin(k2/(n_clad*k0));

% Propagate normally first travel with length of the height of the core and
% cladding layer wrt y offset of the measurement line from comsol.
y       = linspace(0, h*1e-6-y_shift, 400);
Ek_pro  = exp(-1i*ky_clad*y).*(Ek_fil2*ones(size(y)));

% apply Fresnel coefficient for the cladding and air boundary
Ek_pro = apply_Fresnel(Ek_pro, phase_clad, phase_air, n_clad, n_air, 's');

% propagation in air
y2      = linspace(0, d*1e-6, 400);
Ek_pro2 = exp(-1i*ky*y2).*(Ek_pro(:,end)*ones(size(y2)));

Ek_pro  = [Ek_pro Ek_pro2];
Ez_pro  = fftshift(  ifft( fftshift(Ek_pro,1) ),1  );
Ez_tar = Ez_pro(:,end);

% propagate a bit more just for better figure.
y3      = linspace(0, .25*d*1e-6, 400);
Ek_pro3 = exp(-1i*ky*y3).*(Ek_pro(:,end)*ones(size(y3)));

Ek_pro  = [Ek_pro Ek_pro3];
Ez_pro  = fftshift(  ifft( fftshift(Ek_pro,1) ),1  );

%% Checking if original matches simulated
figure(13)
plot(   x2*1e6,     abs(Ez_tar/max(Ez_tar)),    'r',    ...
        x_ori*1e6,  abs(Ez_ori/max(Ez_ori)),    '--b',  ...
        x*1e6,      abs(Ez/max(Ez)),            '-.k')
xlabel('x/ {\mu}m')
ylabel('Normalised |E|')
ylim([0 1.1])
xlim([-60 60])
legend('Simulated', 'Desired Shape', 'Above Grating', 'Location', 'northeast')

%%
writematrix(x2*1e6,                     'fig5b_sim_x.csv');
writematrix(abs(Ez_tar/max(Ez_tar)),    'fig5b_sim_y.csv');

writematrix(x_ori*1e6,                  'fig5b_des_x.csv');
writematrix(abs(Ez_ori/max(Ez_ori)),    'fig5b_des_y.csv');

writematrix(x*1e6,                      'fig5b_gra_x.csv');
writematrix(abs(Ez/max(Ez)),            'fig5b_gra_y.csv');

%%
figure(14)
pcolor(x2*1e6, [y+y_shift y2+h*1e-6 y3+(h+d)*1e-6]*1e6, abs(Ez_pro)');
shading interp

hold on
pcolor(x0, y0, abs(Ez_wg));
shading interp
hold off

xlim([-60 60])
ylim([-h_clad*1e6 inf])
a = colorbar;

xlabel('x/ {\mu}m')
ylabel('y/ {\mu}m')
ylabel(a, '|E|/ kVm^{-1}')

yline(d + h,        '--r', 'LineWidth', 2, 'HandleVisibility','off')
yline(h,            'k', 'LineWidth', 2, 'HandleVisibility','off')
yline(.5*h_core*1e6,'k','LineWidth', 2, 'HandleVisibility','off')

%%
y_new   = [y+y_shift y2+h*1e-6 y3+(h+d)*1e-6]*1e6;
y_max   = y_new(end) - y_new(1) + y0(end) - y0(1);

figure(16)
ax1  = axes;
pl1  = pcolor(x2*1e6, y_new, abs(Ez_pro)');
shading interp

xlim(ax1, [-60 60])
% ylim(ax1, [-h_clad*1e6 inf])

ax2  = axes;
pl2  = pcolor(ax2, x0, y0, abs(Ez_wg));
shading interp
colormap(ax2, 'jet')

ax1.XTick = [];

% Define the size of the plot
ax1_height = .8*(y_new(end) - y_new(1)) / y_max;
ax2_height = .8*(y0(end) - y0(1)) / y_max;

set(ax1, 'Position',[.13 .28 .685 ax1_height]);
set(ax2, 'Position',[.13 .12 .685 ax2_height]);
xlim(ax2, [-60 60])
ylim(ax2, [y0(1) y0(end)])

cb1 = colorbar(ax1,'Position',[.83 .28 .04 ax1_height]);
cb2 = colorbar(ax2,'Position',[.83 .12 .04 ax2_height]);

xlabel(ax2, 'x/ {\mu}m')
ylabel(ax1, 'y/ {\mu}m      ')
ylabel(cb1, '|E|/ kVm^{-1}')

ax1.FontSize = 10;
ax2.FontSize = 10;

yline(ax1, d + h,        '--r', 'LineWidth', 2, 'HandleVisibility','off')
yline(ax1, h,            'k', 'LineWidth', 2, 'HandleVisibility','off')
% yline(ax2, 14,            '-.g', 'LineWidth', 2, 'HandleVisibility','off')
% yline(ax2, .5*h_core*1e6,'k','LineWidth', 2, 'HandleVisibility','off')
%%
writematrix(x2*1e6,         'fig5a_prop_x.csv');
writematrix(y_new,          'fig5a_prop_y.csv');
writematrix(abs(Ez_pro)',   'fig5a_prop_z.csv');

writematrix(x0,             'fig5a_grat_x.csv');
writematrix(y0,             'fig5a_grat_y.csv');
writematrix(abs(Ez_wg),     'fig5a_grat_z.csv');