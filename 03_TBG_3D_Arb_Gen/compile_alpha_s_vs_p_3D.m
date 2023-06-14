% Design of a tilted grating for a 3D arbitrarily shaped output beam with
% controllable polarisation.
%
% Assumptions:
% 1.     target field scalar, from scattering formula match |E| in plane 
% (without the vertical component) to target
% 2.     grating tilt angle is the same everywhere (not correct? local 
% change of grating leads to change of effective tilt angle?)
% 3.     neglect diffraction of pump beam
%
% Design ideas:
% 1.    Conversing two input beams into one square forming a hashtag shape
% grating profile. This might cause some issues with other input pumps
% interact with other respective grating profiles.
% 2.    Having two layers s.t. two input pump doesn't interact with each
% other's grating designs.

%% Parameters
% wavelength
lam     = 1.55e-6; 780e-9;
k0      = 2*pi/lam;
c       = 3e8;

% effective index of fundamental mode
n_eff   = 1.4635;
% propagation constant of mode
beta    = k0*n_eff;
% waist of fundamental mode (vertical)
w0      = 2e-6;
% waist of refractive index profile (vertical)
sig     = 2e-6;
% grating index contrast
dn_g    = 1;            

%% Define target beam (scalar), propagate to grating plane
% grid in grating plane: horizontal cross section at z=0
Lx = 40e-6;     dx = 0.05e-6;
Ly = 40e-6;     dy = 0.05e-6;
x = (-.5*Lx: dx :.5*Lx);
y = (-.5*Ly: dy :.5*Ly);

[xx, yy] = meshgrid(x, y);

n_out = [0, 0, 1];
theta_inc   = linspace(-90,90,2^6)*pi/180;
alpha_s     = zeros(1, length(theta_inc));
alpha_p     = alpha_s;

E1 = define_target_field(xx, yy, beta, n_out,...
    'plane', 2.5e-6, 50e-6);
[Lam_grat0_1, alp_grat0_1, alp_tilt0_1] = compute_grating_angles(...
    n_out(1), n_out(2), n_out(3), [1, 0, 0]);
Lam_grat0_1 = Lam_grat0_1 * lam/n_eff;

% =1 horizontal input polarisation; =2 vertical input polarisation
pol = 1;
% transverse prop. constant of pump
ky = 0;     

for ite2 = 1:length(theta_inc)
    [al_s, al_p] = compute_alpha(...
        Lam_grat0_1, theta_inc(ite2), alp_tilt0_1,...
        w0, sig, beta, dn_g, n_eff, ky);
    alpha_s(ite2) = al_s;
    alpha_p(ite2) = al_p;
end

%% plot of scattering rate vs K_pump(theta_inc)
figure(2); clf;
plot(...
    theta_inc*180/pi, alpha_s,...
    theta_inc*180/pi, alpha_p);
xlabel('rotated angle, \theta_{inc} / [deg]'); ylabel('scattering rate');
legend('\alpha_s','\alpha_p');