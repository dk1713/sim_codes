%% Design and test tilted Gaussian beam

w0      = 3e-6;             % beam waist
phi_deg = 70;               % tilt angle of Gaussian in deg
phi     = phi_deg * pi/180; % tilt angle of Gaussian (90deg = vertical)
angle   = pi - phi;

d       = 100e-6;           % focal distance
lam     = .78e-6;           % wavelength

n_clad  = 1.447;
k0      = 2*pi*n_clad/lam;

y0      = pi * w0^2 * n_clad / lam;    % Rayleigh range

% on the focus centred at the focus
x       = linspace(-60e-6, 60e-6, 5000)';
x_rot   = x ./cos(.5*pi - phi);
y_rot   = - x .*tan(.5*pi - phi);


wz      = w0*sqrt(1+(y_rot/y0).^2);
Rz      = 1./y_rot .* (y_rot.^2+y0^2);
etaz    = 0.5*atan(y_rot/y0);

E       = sqrt(w0./wz) .* exp(-1i*k0*y_rot + 1i*etaz ...
                        - x_rot.^2./wz.^2 ...
                        - 1i*k0*x_rot.^2./(2*Rz));
Ek      = fftshift( fft( fftshift(E) ) );

dat = [60+1e6*x' E.'];

figure(1)
plot(x, abs(E))
xlabel('x')
ylabel('y')

%% Tilted angle
k       = linspace(-pi, pi, length(x) )'./(x(1) - x(2));
ky      = real( sqrt(k0^2 - k.^2) );

% find angle of beam propagation
ii      = find(ky>0);
Sx      = sum(abs(Ek(ii)).^2 .* k(ii));
Sy      = sum(abs(Ek(ii)).^2 .* ky(ii));

phi     = pi - atan2(Sy, Sx);

x_focus = -d * cos(phi);    % x-coor of focus
y_focus = d * sin(phi);     % y-coor of focus

disp('Gaussian beam design:')
disp(['   Waist (um)            : ' num2str(1e6*w0)])
disp(['   Waist position (um)   : ' num2str(1e6*d)])
disp(['   Rayleigh range (um)   : ' num2str(1e6*y0)])
disp(['   Tilt angle (deg)      : ' num2str(phi*180/pi)])

%% Propagate the wave
y       = linspace(-100, 100, 500)*1e-6;
Ek_pro  = exp(-1i*ky*y).*(Ek*ones(size(y)));
Ez_pro  = fftshift(  ifft( fftshift(Ek_pro,1) ),1  );

figure(2)
pcolor(x, y, abs(Ez_pro)')
xlabel('x')
ylabel('y')
shading flat
colorbar
hold on
xline(0, '--r')
yline(0, '--r')
hold off
axis equal

figure(3)
pcolor(x + x_focus, y + y_focus, abs(Ez_pro)')
xlabel('x')
ylabel('y')
shading flat
colorbar
hold on
xline(0, '--r')
yline(0, '--r')
hold off
axis equal

%% E on the grating
Ek_grat = exp(-1i * ky * -y_focus).*(Ek);
Ez_grat = fftshift(  ifft( fftshift(Ek_grat) ) );

figure(4)
plot(x + x_focus, abs(Ez_grat))
xlabel('x')
ylabel('y')

%% Parameters for grating
% exp term
Ez_amp      = abs(Ez_grat);
Ez_exp_img  = imag(log(Ez_grat));

figure(11)
plot(x, Ez_exp_img)
xlabel('x')
ylabel('Img part of exp comp')
yline(pi,'--r')
yline(-pi, '--r')