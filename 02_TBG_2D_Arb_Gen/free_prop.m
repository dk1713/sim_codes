%% free space propagation and Gaussian beam fitting
% load electric field above waveguide from Comsol simulation
% Parameters
lam     = 0.78e-6;
n_clad  = 1.447;
dn      = 5e-3;

% Target distance and beam waist
d       = 100;
w_0     = 15;
phi_deg = 80;
n       = 8;
outP    = 0.05;

load([  'data/bragg_grating',     ...
        '_dimensions_',         num2str(20), 'x', num2str(130), ...
        '_grat_length_',        num2str(120),                   ...
        '_lambda_',             num2str(lam * 1e6),              ...
        '_n_clad_',             num2str(n_clad),                ...
        '_dn_',                 num2str(dn),                    ...
        '_scattered_tilt_',     num2str(phi_deg),               ...
        '_target_dist_',        num2str(d),                     ...
        '_target_waist_',       num2str(w_0),                   ...
        '_gausssian_order_',    num2str(n),                     ...
        '_power_ratio_',        num2str(outP),                  ...
        '.mat']);

%% Initialisation
x       = x * 1e-6;
k0      = 2*pi * n_clad / lam;
y_shift = 9.5e-6;

% Defining the element size
N       = length(x);
dx      = x(2) - x(1);

% Defining the Fourice space
k       = linspace(-pi,pi,N)'/dx;
Ek      = fftshift( fft( fftshift(Ez) ) );

%% Filter the spectrum to get the side-scattered light only
lower_bound = -.8*k0;
upper_bound = k0;

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
figure(1)
plot(   x,  abs(Ez),        'b', ...
        x2, abs(Ez_fil2),   '--r' )
xlabel('x')
ylabel('|E|')
legend('original', 'filtered')

figure(2)
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

% Propagate normally
y       = linspace(-0,200,500)*1e-6;
Ek_pro  = exp(-1i*ky*y).*(Ek_fil2*ones(size(y)));
Ez_pro  = fftshift(  ifft( fftshift(Ek_pro,1) ),1  );

figure(3)
pcolor(x2, y + y_shift, abs(Ez_pro)')
xlabel('x')
ylabel('y')
shading flat
colorbar
axis equal
%%
y_tar   = d*1e-6 - y_shift;
Ek_tar  = exp(-1i*ky*y_tar).*Ek_fil2;
Ez_tar  = fftshift(ifft( fftshift(Ek_tar)) );

figure(4)
plot(   x_coor, abs(E_tar/max(E_tar)), 'b',    ...
        x2, abs(Ez_tar/max(Ez_tar)),    '--r' )
xlabel('x')
ylabel('y')
ylim([0 1.1])
legend('Original', 'Propagated')
title('field on target')