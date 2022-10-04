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
n       = 1;
k0      = 2*pi*n_air/lambda;
%% Target specification
% distance from the top of the chip.
dist    = 50e-6;
% estimated waist at the surface
waist   = 5e-6;
%% On the Surface
% 3D Gaussian on the Surface
% init:
psi     = 30*pi/180;

L_x         = 2*w_x;
L_y         = 2*w_y;

Nx          = 2^10; 
Ny          = Nx/2;
dx          = L_x / Nx;
dy          = L_y / Ny;

EE      = zeros(Ny,Nx);
kk_z    = zeros(Ny,Nx);