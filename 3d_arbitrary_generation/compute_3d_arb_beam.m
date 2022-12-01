%% Computing grating profile for 3D Arbitrary beam generation
% Code to create tilted Bragg grating profile for 3D arbitrary beam
% generation. 2 steps are needed and this codes is 2nd part out of 2 steps.
%       1. Define the required field shape on the surface or propagate the
%       required field shape onto the surface.
%       2. Propagate the field on the surface onto the grating plane.
%       3. Loop through the vertical layers parallel to the slab mode
%       coming from the initial grating. This considers the model where the
%       dispersion is not assumed.
%
%% Device specification
% Indices [1]
n_air   = 1;
n_clad  = 1.4555;
n_core  = 1.4608;
n_eff   = 1.4635;
dn_g    = 3e-3;

% Heights of the layers [m]
h_clad  = 15e-6;
h_core  = 5e-6;

lam     = 780e-9;
k0      = 2*pi*n_air/lam;
%% Target specification
% distance from the top of the chip.
pos_tar     = [-30e-6, 20e-6, 100e-6];

% compute for angles in grating profiles:
[theta_grat, theta_tilt] = grating_angles(pos_tar, k0);

% estimated waist at the surface
waist_tar   = [5e-6, 10e-6];
%% Init
% 1. defining the domain
L_x     = 200e-6;
L_y     = 200e-6;

N_x     = 2^10; 
N_y     = N_x/2^6;

x       = linspace(-.5*L_x, .5*L_x, N_x);
y       = linspace(-.5*L_y, .5*L_y, N_y);

% 2. defining the arbitrary beam
n       = 1;
[xx, yy]    = meshgrid(x, y);
k_cen_x     = k0*pos_tar(1)/sqrt(pos_tar(1)^2+pos_tar(2)^2+pos_tar(3)^2);
k_cen_y     = k0*pos_tar(2)/sqrt(pos_tar(1)^2+pos_tar(2)^2+pos_tar(3)^2);

EE        = exp( -(...
    (.5*(xx - pos_tar(1))./waist_tar(1)).^2     ...
    + (.5*(yy - pos_tar(2))./waist_tar(2)).^2   ...
    ).^n );

figure(1)
pcolor(xx, yy, abs(EE).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal
xline(pos_tar(1), 'r');
yline(pos_tar(2), 'r');

%% Fourier space
dk_x    = 2*pi/L_x;
dk_y    = 2*pi/L_y;
k_x     = (-N_x/2:N_x/2-1) * dk_x;
k_y     = (-N_y/2:N_y/2-1) * dk_y;
[kk_x, kk_y] = meshgrid(k_x, k_y);

% computing for wave vectors for air
kk_t    = sqrt((kk_x + k_cen_x).^2 + (kk_y + k_cen_y).^2);
kk_z    = real( sqrt(k0^2 - kk_t.^2) );

% checking for limits of central freqency
fprintf('--------------------------------------------------------------\n')
fprintf('Checking for limits of central wavenumber\n')
fprintf('k_t(x) range: (%3.3e, %3.3e)',min(k_x+k_cen_x),max(k_x+k_cen_x))
fprintf('\nk_cen_x  = %3.3e\n', k_cen_x)
fprintf('k_t(y) range: (%3.3e, %3.3e)',min(k_y+k_cen_y),max(k_y+k_cen_y))
fprintf('\nk_cen_y  = %3.3e\n', k_cen_y)
fprintf('--------------------------------------------------------------\n')

EE_k    = fftshift(fft2(fftshift(EE))) ...
    .* exp(1i*kk_z.*-1*pos_tar(3)); % propagation back onto the surface.
EE      = fftshift(ifft2(fftshift(EE_k)))...
    .* exp(1i*(k_cen_x.*xx + k_cen_y.*yy)); % add back central frequency.

figure(2)
pcolor(xx, yy, abs(EE).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal
xline(pos_tar(1), 'r');yline(pos_tar(2), 'r');
xline(0, 'k');yline(0, 'k');
xline(L_x/4,'g');yline(L_y/4,'g');xline(-L_x/4,'g');yline(-L_y/4,'g');

%% Select the domain of the grating surface (Zooming in)
Lg_x    = 100e-6;
Lg_y    = 100e-6;

Ng_x    = N_x/ceil(L_x/Lg_x);
Ng_y    = N_y/ceil(L_y/Lg_y);

index_x = (0:Ng_x-1) + find(x>-Lg_x/2,1);
index_y = (0:Ng_y-1) + find(y>-Lg_y/2,1);

xg      = x(index_x);
yg      = y(index_y);

[xx_g, yy_g] = meshgrid(xg, yg);
EE_g    = EE(index_y,index_x);

figure(3)
pcolor(xx_g, yy_g, abs(EE_g).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal

%% Propagate the field onto the grating plane
% re-evaluating the Fourier domain
dk_g_x  = 2*pi/Lg_x;
dk_g_y  = 2*pi/Lg_y;
k_g_x   = (-Ng_x/2:Ng_x/2-1) * dk_g_x;
k_g_y   = (-Ng_y/2:Ng_y/2-1) * dk_g_y;
[kk_g_x, kk_g_y] = meshgrid(k_g_x, k_g_y);

% The grating domain is fixed hence xx_g and yy_g is carried over.
% Propagating from surface to 1st boundary between cladding and core.
[EE_g, k_cen_x, k_cen_y] = prop_boundary(EE_g, xx_g, yy_g, ...
    k0*n_air, k0*n_clad, k_cen_x, k_cen_y, ...
    kk_g_x, kk_g_y, n_air, n_clad, -h_clad, 's');

% new k_cen's are calculated each propagation.

% Propagating from 1st boundary to grating plane.
[EE_g, k_cen_x, k_cen_y] = prop_boundary(EE_g, xx_g, yy_g, ...
    k0*n_clad, k0*n_core, k_cen_x, k_cen_y, ...
    kk_g_x, kk_g_y, n_clad, n_core, -h_core, 's');

figure(4)
pcolor(xx_g, yy_g, abs(EE_g).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal

%% Compute for power ratios
% Need to consider
% 1. Integration over length, x
% 2. diffraction angle

% init
k_core  = k0*n_core;
k_t     = sqrt( k_cen_x.^2 + k_cen_y.^2 ); % need to update to make it general later by finding the 1st derivative!
k_z     = real( sqrt(k_core^2 - k_t.^2) );
% computing for diffraction angles (propagation angle in core layer)
phi     = pi/2 - sign(k_cen_x)*acos(k_z/k_core);

power_geo   = abs(EE_g.^2) .* sin(phi);
power_ratio = trapz(xx_g(1,:), power_geo, 2);
power_ratio = power_ratio/max(power_ratio);

figure(5)
plot(yy_g(:,1), power_ratio, 'x', 'markersize', 10);

%% Compute for constant of normalisation
% init memo allo
dn_gs   = zeros(size(EE_g));
period  = zeros(size(EE_g));
efficiency = zeros(size(power_ratio));
beta    = k0*n_eff;

% parameters in BTA
w_0     = 2e-6; % may need to change.
sigma   = .5*h_core;
w       = w_0*sigma/sqrt(w_0^2 + sigma^2);

% computing for incident angle relative to grating slices.
k_out   = [k_cen_x, k_cen_y, k_z];
k_in    = [beta, 0, 0];
theta_inc = .5*acos(dot(k_out,-k_in)/norm(k_out)/norm(k_in));

figure(10)
fprintf('computing for optimum dn_g < %2.4e... \n', dn_g);
for i = 1:length(power_ratio)
    E_grat  = EE_g(i,:)*power_ratio(i);
    
    phase   = unwrap(angle(E_grat));
    dphase  = beta + c_diff(xg, phase);
    Lam     = 2*pi./abs(dphase);
    period(i,:) = Lam;
    Pz_amp  = abs(E_grat).^2;

    % init for the loop
    F       = griddedInterpolant(xg, Pz_amp, 'spline');
    fun     = @(x) F(x);
    C       = 1/integral(fun, min(xg), max(xg));
    max_dng = 1e2;
    eta     = 1;
    
    denu = zeros(size(xg));
    for ii = 1:length(xg)
        denu(ii) = abs(1 - eta*C * integral(fun, xg(1), xg(ii))); % may need abs to improve accuracy in error function.
    end

    kpump_n     = k_cen_x*cos(theta_grat) + k_cen_y*sin(theta_grat);
    dng_amp     = sqrt(eta*C * fun(xg) ./ denu);
    dng_amp     = 2*cos(theta_inc).^2 .* sin(2*theta_tilt) .* dng_amp;
    dng_amp     = n_eff .* Lam ./ w .* sqrt(w_0 ./ sqrt(pi)) .* dng_amp;

    dng_amp     = dng_amp ./ exp(-1*(w./sin(2*theta_tilt)).^2 ...
        .*(kpump_n.*cos(theta_tilt).^2-pi./Lam).^2);

    dng_temp    = dng_amp ./ (1 + cos(2*theta_inc).^2) ./sqrt(sin(phi));
    
    hold on
    plot(xg, dng_temp);
%     plot(xg, 1./exp(-1*(w./sin(2*theta_tilt)).^2 ...
%         .*(kpump_n.*cos(theta_tilt).^2-pi./Lam).^2));
    xlabel('x');
    ylabel('intensity');

%     while 1e3*abs(dn_g - max_dng) > 1e-4
%         if eta < 0
%             fprintf('eta = %2.2f has to be positive value!\n', eta);
%             break
%         end
%         
%         denu = zeros(size(xg));
%         for ii = 1:length(xg)
%             denu(ii) = abs(1 - eta*C * integral(fun, 10*xg(1), xg(ii))); % may need abs to improve accuracy in error function.
%         end
% 
%         kpump_n     = k_cen_x*cos(theta_grat) + k_cen_y*sin(theta_grat);
%         dng_amp     = sqrt(eta*C * fun(xg) ./ denu);
%         dng_amp     = 2*cos(theta_inc).^2 .* sin(2*theta_tilt) .* dng_amp;
%         dng_amp     = n_eff .* Lam ./ w .* sqrt(w_0 ./ sqrt(pi)) .* dng_amp;
%         
%         disp(max(dng_amp));
%         
%         dng_amp     = dng_amp ./ exp(-2*(w./sin(2*theta_tilt)).^2 ...
%             .*(kpump_n.*cos(theta_tilt).^2-pi./Lam).^2);
%         
%         disp(max(exp(-2*(w./sin(2*theta_tilt)).^2 ...
%             .*(kpump_n.*cos(theta_tilt).^2-pi./Lam).^2))); %%TODO problem area
%         
%         dng_temp    = dng_amp ./ (1 + cos(2*theta_inc).^2)./sqrt(sin(phi));
%         
%         dn_gs(i,:)  = dng_temp;
%         max_dng = max(dng_temp);
%         dng_diff = 1e3*(dn_g - max_dng);
%     break;
%         if dng_diff < 0
%             eta = eta - .05;
%         elseif dng_diff > 1
%             eta = eta + .05;
%         else
%             eta = eta + .01*dng_diff;
%         end
%     end
%     break;
%     alpha_ana = sqrt(2*pi)*pi^2*w.^2*dng_temp.^2 ./ ...
%         w_0 ./ (2*cos(theta_inc)^2)^2 ./ Lam ./ sin(theta_tilt).^2 ...
%         ./ n_eff^2 .* sin(phi).* ...
%         exp(-2*(w./sin(2*theta_tilt)).^2 ...
%         .*(kpump_n.*cos(theta_tilt).^2-pi./Lam).^2).^2;
%     efficiency(i) = 100 - 100*exp(-trapz(xg, alpha_ana));
end
hold off

% saving data needed:
EE_grat = exp(-1i*beta*xg) .* EE_g;
phase   = angle(EE_grat);
% Required .mat files
save('3d_gauss.mat', 'xg', 'yg', 'EE_grat', 'phase', 'period', 'dn_gs', 'efficiency');