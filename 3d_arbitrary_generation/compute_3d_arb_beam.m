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
dn_g    = 1e-3;

% Heights of the layers [m]
h_clad  = 15e-6;
h_core  = 5e-6;

lam     = 780e-9;
k0      = 2*pi*n_air/lam;
%% Target specification
% distance from the top of the chip.
pos_tar     = [20e-6, 10e-6, 50e-6];
% estimated waist at the surface
waist_tar   = [5e-6, 10e-6];
%% 3D Gaussian on the Surface
% 1. defining the domain
L_x     = 200e-6;
L_y     = 200e-6;

N_x     = 2^12; 
N_y     = N_x/(2^7);

x       = linspace(-.5*L_x, .5*L_x, N_x);
y       = linspace(-.5*L_y, .5*L_y, N_y);

% 2. defining the arbitrary beam
n       = 1;
[xx, yy]    = meshgrid(x, y);
EE_1        = exp( -(...
    (.5*(xx - pos_tar(1))./waist_tar(1)).^2     ...
    + (.5*(yy - pos_tar(2))./waist_tar(2)).^2   ...
    ).^n );

figure(1)
pcolor(xx, yy, abs(EE_1).^2)
shading flat
xlabel('x')
ylabel('y')
title('Intensity on surface')
colorbar
axis equal