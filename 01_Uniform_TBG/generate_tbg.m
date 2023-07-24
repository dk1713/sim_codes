%% Generation of the dataset from COMSOL simulations
% Preset COMSOL file: 'tbg_constant_dng.mph'
model = mphload('tbg_constant_dng.mph');
% This script generates the numerical dataset from the COMSOL simulation.
% Note that this takes time and might need to check if preset COMSOL file
% works beforehand.

%% Inside the COMSOL file          
%    -------------power_up--------------
%    |                                 |
% power_in                          power_out
%    |                                 |
%    ------------power_down-------------
% The four data points are taken. Power coming in (input mode) and power
% leakage and power flowed out. From tilted Bragg grating, we mainly want
% power_up and power_out.
%
% power_in is set as 1W in COMSOL and power up and down is only taken when
% tilted Bragg grating scatters light at 90[deg] (theta = 45[deg]).

%% Taking up parameter data from .mph
% for labelling
L       = mphglobal(model, 'L');
H       = mphglobal(model, 'H');
L_g     = mphglobal(model, 'L_g');
H_g     = mphglobal(model, 'H_g');
n_cl    = mphglobal(model, 'n_cl');
n_co    = mphglobal(model, 'n_co');
sigma   = mphglobal(model, 'sigma');
dn      = mphglobal(model, 'dn');
% Note that the dimensions are in um so careful with the unit.

%% init
dn_g    = 4e-3;
lam     = 780;
phi_deg = 120;
% Default values:
model.param('par3').set( 'dn_g',    num2str(dn_g) );
model.param('par2').set( 'lambda0', [num2str(lam), ' [nm]']);
model.param('par2').set( 'phi',     [num2str(phi_deg), ' [deg]'] );

model_spec = [                      ...
    '_L_',      num2str(L),         ...
    '_H_',      num2str(H),         ...
    '_L_g_',    num2str(L_g),       ...
    '_H_g_',    num2str(H_g),       ...
    '_n_cl_',   num2str(n_cl),      ...
    '_n_co_',   num2str(n_co),      ...
    '_sigma_',  num2str(sigma),     ...
    '_dn_',     num2str(dn),        ...
    '_phi_',    num2str(phi_deg),   ...
    ];

%% Input Mode
disp('v---------------------Computation started---------------------v');
disp('Computing for width of the input mode');
% setting no dng s.t. there is no back reflection. 
model.param('par3').set( 'dn_g',   num2str(0) );
model.study('std1').run
n_eff   = mphglobal(model, 'real(ewfd.neff_1)');

y       = linspace(-.5*H, .5*H, 2^9);
x       = -.5*L;
xx      = [x*ones(size(y)); y];
normE   = mphinterp(model, 'ewfd.normE', 'coord', xx);
Poavx   = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
power_in= trapz(y*1e-6, Poavx);

% function handles
gauss       =   @(x, a, b, c)...
    a*exp(-((x-b)/c).^2);
fit_gauss   =   @(p, x, data)...
    sqrt(mean( (gauss(x,p(1),p(2),p(3)) - data).^2  ));

ooo = optimset('TolX',1e-9);
[par,~] = fminsearch(@(p) fit_gauss(p, y, normE), ...
        [max(normE) sum(normE.*y)/sum(normE) 1], ooo);

fprintf('Gaussian width is %2.4f [um]\n', par(3))
w_0     = par(3);
theta   = mphglobal(model, 'theta');

% Checking for input mode
plot_check = 0;
if plot_check
    figure(1); clf;
    plot(y, normE, 'b', y, gauss(y, par(1), par(2), par(3)), '--r')
    xlabel('x [{\mu}m]')
    ylabel('y [V / m]')
    legend('numerical', 'fitting' )
end
disp('^--------------------Computation completed--------------------^');

%% Init for three parametric sweep
total_size  = 2^9;
y       = linspace(-.5*H, .5*H, total_size);
x       = .5*L;
xx      = [x*ones(size(y)); y];

%% dn_g efficiency studies
disp('v---------------------Computation started---------------------v');
fprintf('Computing for index Case: dn = %2.3e and varying dn_g\n', dn);

dn_gs           = linspace(1, 5, 9) * 1e-3;
num_of_ite      = length(dn_gs);
power_out_dng   = ones(size(dn_gs));

for iter = 1:num_of_ite
    fprintf('Computing for dng = %2.3e', dn_gs(iter));
    fprintf('    (iteration #%i of %i iterations)\n',...
        iter, num_of_ite);

    model.param('par3').set('dn_g', num2str(dn_gs(iter)));
    model.study('std1').run;

    % Measuring time-average powerflow in y
    Poavx = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
    power_out_dng(iter) = trapz(y*1e-6, Poavx);
end
disp('^--------------------Computation completed--------------------^');

%% lambdas efficiency studies
disp('v---------------------Computation started---------------------v');
fprintf('Computing for index Case: dn = %2.3e and varying lambda\n', dn);

% reset to default
model.param('par3').set( 'dn_g',   num2str(dn_g) );

lambdas         = linspace(740, 820, 9);
num_of_ite      = length(lambdas);
power_out_lam   = ones(size(lambdas));

for iter = 1:num_of_ite
    fprintf('Computing for lambda = %4.1f [nm]', lambdas(iter));
    fprintf('    (iteration #%i of %i iterations)\n', ...
        iter, num_of_ite);
    
    model.param('par2').set('lambda0', [num2str(lambdas(iter)), ' [nm]']);
    model.study('std1').run;
    
    % Measuring time-average powerflow in y
    Poavx = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
    power_out_lam(iter) = trapz(y*1e-6, Poavx);
end
disp('^--------------------Computation completed--------------------^');

%% tilt angle efficiency studies
disp('v---------------------Computation started---------------------v');
fprintf('Computing for index Case: dn = %2.3e and varying theta\n', dn);

% reset to default
model.param('par2').set( 'lambda0',   [num2str(lam), ' [nm]']);

theta_deg   = .5*phi_deg;
thetas      = linspace(theta_deg - 4, theta_deg + 4, 9);
num_of_ite = length(thetas);
power_out_the  = ones(size(thetas));

for iter = 1:num_of_ite
    fprintf('Computing for theta = %4.1f [deg]', thetas(iter));
    fprintf('    (iteration #%i of %i iterations)\n', ...
        iter, num_of_ite);

    model.param('par2').set('theta', [num2str(thetas(iter)), ' [deg]']);
    model.study('std1').run;
    
    % Measuring time-average powerflow in y
    Poavx = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
    power_out_the(iter) = trapz(y*1e-6, Poavx);
end
disp('^--------------------Computation completed--------------------^');

%% saving the dataset
filename    =   ['data/tbg_efficiency'  model_spec '.mat'];
save(filename, 'n_eff', 'w_0', 'power_in', 'theta',  ...
    'dn_gs',    'power_out_dng',           ...
    'lambdas',  'power_out_lam',        ...
    'thetas',   'power_out_the');