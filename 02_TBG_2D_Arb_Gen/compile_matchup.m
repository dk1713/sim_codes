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
lam     = 780e-9;

% Target distance and beam waist [um]
% r       = 30;
r       = [-30, 30];
% r       = [20, 30] + 10;
d       = 50;
w_0     = 3;
n       = 1;
outP    = 0.2;

h       = (h_clad + .5*h_core)*1e6;

% Dimensions
L       = 120;
H       = 30;
L_g     = 100;

% tar_xs string generator
newStr = num2str(r(1));
for ite = 2:length(r)
    newStr = [newStr 'x' num2str(r(ite))]; %#ok<AGROW>
end

load([  'data/grating_design',                                  ...
        '_dim_',                num2str(H), 'x', num2str(L),    ...
        '_grat_len_',           num2str(L_g),                   ...
        '_lambda_',             num2str(lam * 1e6),             ...
        '_n_clad_',             num2str(n_clad),                ...
        '_dn_',                 num2str(dn),                    ...
        '_tar_dist_',           num2str(d),                     ...
        '_tar_xs_',             newStr,                         ...
        '_tar_waist_',          num2str(w_0),                   ...
        '_gauss_order_',        num2str(n),                     ...
        '_eta_',                num2str(outP),                  ...
        '.mat']);
%% original figure (only for Gaussian)
% figure(10); clf; colororder(["#8040E6";"#1AA640";"#E68000"]);
% plot(x_ori*1e6,  abs(Ez_ori/max(Ez_ori)), '-', 'linewidth', 1)
% xlabel('x/ {\mu}m')
% ylabel('Normalised |E|')
% ylim([0 1.1])
% xlim([-15 15])
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
figure(11); clf; colororder(["#8040E6";"#1AA640";"#E68000"]);
plot(   x*1e6,  abs(Ez),         ...
        x2*1e6, abs(Ez_fil2),   '--', 'linewidth', 1)
xlabel('x / {\mu}m')
ylabel('|E| / kVm^{-1}')
xline(-50,  '--k','HandleVisibility','off')
xline(50,  '--k','HandleVisibility','off')
legend('original', 'filtered')

figure(12); clf; colororder(["#8040E6";"#1AA640";"#E68000"]);
plot(   k,  abs(Ek),         ...
        k2, abs(Ek_fil2),   '--', 'linewidth', 1)
xlabel('k')
ylabel('|Ek|')
legend('original', 'filtered')
hold on
xline(lower_bound,  '--k','HandleVisibility','off')
xline(upper_bound,  '--k','HandleVisibility','off')
hold off
xlim(3*[lower_bound,upper_bound])

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
figure(13); clf; colororder(["#8040E6";"#1AA640";"#E68000"]);
plot(   x2*1e6,     abs(Ez_tar/max(Ez_tar)),     ...
        x_ori*1e6,  abs(Ez_ori/max(Ez_ori)),    '--',  ...
        x*1e6,      abs(Ez/max(Ez)),            '-.', 'linewidth', 1)
xlabel('x/ {\mu}m')
ylabel('Normalised |E|')
ylim([0 1.1])
xlim([-60 60])
legend('Simulated', 'Desired Shape', 'Above Grating', 'Location', 'northwest')

%% plot of index distribution
figure(14); clf;
colororder(["#8040E6";"#1AA640"])
yyaxis left
plot( x_ori*1e6, dng_need*1e3, 'linewidth', 1);
xlabel('x/ {\mu}m');
ylabel('index modulation, {\Delta}n_g / [10^{-3}]')

yyaxis right
plot( x_ori*1e6, phase, 'linewidth', 1 );
ylabel('required phase, {\Phi}(x)')
xline(-50,  '--k','HandleVisibility','off')
xline(50,  '--k','HandleVisibility','off')

%%
% figure(15); clf;
% pcolor(x2*1e6, [y+y_shift y2+h*1e-6 y3+(h+d)*1e-6]*1e6, abs(Ez_pro)');
% shading interp
% 
% hold on
% pcolor(x0, y0, abs(Ez_wg));
% shading interp
% hold off
% 
% xlim([-60 60])
% ylim([-h_clad*1e6 inf])
% a = colorbar;
% 
% xlabel('x/ {\mu}m');
% ylabel('y/ {\mu}m');
% ylabel(a, '|E|/ kVm^{-1}');
% 
% yline(d + h,        '--r', 'LineWidth', 2, 'HandleVisibility','off')
% yline(.5*H,            'k', 'LineWidth', 2, 'HandleVisibility','off')
% % yline(.5*h_core*1e6,'--k','LineWidth', 1, 'HandleVisibility','off')

%%
y_new   = [y+y_shift y2+h*1e-6 y3+(h+d)*1e-6]*1e6;
y_max   = y_new(end) - y_new(1) + y0(end) - y0(1);

figure(16); clf;
ax1  = axes;
pl1  = pcolor(x2*1e6, y_new, abs(Ez_pro)');
shading interp

xlim(ax1, [-60 60])
% ylim(ax1, [-h_clad*1e6 inf])

ax2  = axes;
pl2  = pcolor(ax2, x0, y0, abs(Ez_wg));
shading interp
% colormap(ax2, 'jet')

ax1.XTick = [];

% Define the size of the plot
ax1_height = .8*(y_new(end) - y_new(1)) / y_max;
ax2_height = .8*(y0(end) - y0(1)) / y_max;

set(ax1, 'Position',[.13 .38 .685 ax1_height]);
set(ax2, 'Position',[.13 .12 .685 ax2_height]);
xlim(ax2, [-60 60])
ylim(ax2, [y0(1) y0(end)])

cb1 = colorbar(ax1,'Position',[.83 .38 .04 ax1_height]);
cb2 = colorbar(ax2,'Position',[.83 .12 .04 ax2_height]);

xlabel(ax2, 'x/ {\mu}m')
ylabel(ax1, 'y/ {\mu}m      ')
ylabel(cb1, '|E|/ kVm^{-1}')

ax1.FontSize = 10;
ax2.FontSize = 10;

yline(ax1, d + h,        '--r', 'LineWidth', 2, 'HandleVisibility','off')
% yline(ax1, h,            'k', 'LineWidth', 2, 'HandleVisibility','off')
% yline(ax2, 14,            '-.g', 'LineWidth', 2, 'HandleVisibility','off')
% yline(ax2, .5*h_core*1e6,'k','LineWidth', 2, 'HandleVisibility','off')

%% plot of the close-up of the grating feature
figure(17); clf;
pcolor(x0, y0, abs(Ez_wg));
shading interp
% colormap(ax2, 'jet')
ylim([-5, 5]);
xlim([-15, 15]);
xlabel('x/ {\mu}m')
ylabel('y/ {\mu}m');
cb1 = colorbar;
ylabel(cb1, '|E|/ kVm^{-1}')
axis equal