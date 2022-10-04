%% Efficiency of an ion trap by beam tracing method
% Compiling the final figures from the data generated.

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
%% Target specification
% distance from the top of the chip.
dist    = 8e-3;
% estimated waist at the surface
waist   = 2e-3;
% x position of the target
xs_tar = [100e-6 -100e-6];

%% Filename
model_spec = [                              ...
    '_n_clad_',     n_clad,                 ...
    '_n_core_',     n_core,                 ...
    '_h_clad_',     h_clad,                 ...
    '_h_core_',     h_core,                 ...
    '_lam_',        lam,                    ...
    '_dist_',       dist,                   ...
    '_waist_',      waist,                  ...
    '_spacing_',    xs_tar(2) - xs_tar(1),  ...
    '_num_g_',      length(xs_tar),         ...
    '_max_dng_',    dn_g,                   ...
    ];

filename    =   ['data/tbg_multi_gauss'  model_spec '.mat'];
load(filename)

%% final figures
figure(1)
colororder({'b','r'})
yyaxis left
plot(1e3*x, 1e9*period, 'b')
yline(1e9*lam/(1*cos(phi+.5*pi) + n_eff), '--k', 'period for constant k');
xlabel('x / [mm]')
ylabel('grating period, {\Lambda} / [nm]')

yyaxis right
plot(1e3*x, dn_gs, 'r');
xlabel('x / [mm]')
ylabel('index modulation, {\Delta}n_g')