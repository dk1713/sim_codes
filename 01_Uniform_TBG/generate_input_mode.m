%% Generation of the input mode analysis dataset from COMSOL simulations
% Preset COMSOL file: 'tbg_constant_dng.mph'
model = mphload('tbg_constant_dng.mph');
% This script generates the numerical dataset from the COMSOL simulation.
% Note that this takes time and might need to check if preset COMSOL file
% works beforehand.
%
% It uses the same COMSOL file for uniform grating analysis to study the
% input mode is affected by the sigma and also looking into the
% efficiencies when scattering upwards at 90 degs.

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
dn_g    = 3e-3;
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
normE0  = mphinterp(model, 'ewfd.normE', 'coord', xx);

% function handles
gauss       =   @(x, a, b, c)...
    a*exp(-((x-b)/c).^2);
fit_gauss   =   @(p, x, data)...
    sqrt(mean( (gauss(x,p(1),p(2),p(3)) - data).^2  ));

ooo = optimset('TolX',1e-9);
[par0,~] = fminsearch(@(p) fit_gauss(p, y, normE0), ...
        [max(normE0) sum(normE0.*y)/sum(normE0) 1], ooo);

fprintf('Gaussian width is %2.4f [um]\n', par0(3))
w_0     = par0(3);
theta   = mphglobal(model, 'theta');
fitting0= gauss(y, par0(1), par0(2), par0(3));

% Checking for input mode
figure(1); clf;
plot(y, normE0, '-', y, fitting0, '--');
xlabel('x [{\mu}m]');
ylabel('y [V / m]');
legend('numerical', 'fitting' );

error0 = sum((gauss(y,par0(1),par0(2),par0(3))-normE0).^2) ...
    /sum((normE0 - mean(normE0)).^2);
fprintf('Relative squared error = %2.4e \n', error0);

disp('^--------------------Computation completed--------------------^');

%% Parametric sweep for sigma
disp('v---------------------Computation started---------------------v');
fprintf('Computing for n_eff and w_0 for varying sigma\n');

% reset to default
model.param('par3').set( 'dn_g',   num2str(dn_g) );

sigmas          = linspace(0.5, 5, 10);
num_of_ite      = length(sigmas);
w_0s            = ones(size(sigmas));
n_effs          = ones(size(sigmas));
power_out_sig   = ones(size(sigmas));

errors          = ones(size(sigmas));
normEs          = ones(length(sigmas), length(y));
fittings        = ones(length(sigmas), length(y));
% errors          = ones(1,2);
% normEs          = ones(2,length(y));
% fittings        = ones(2,length(y));
% err_ite         = 0;

for iter = 1:num_of_ite
    fprintf('Computing for sigma = %4.1f [um]', sigmas(iter));
    fprintf('    (iteration #%i of %i iterations)\n', ...
        iter, num_of_ite);
    
    model.param('par3').set('sigma', [num2str(sigmas(iter)), ' [um]']);
    model.study('std1').run;
    
    % output
    x       = .5*L;
    xx      = [x*ones(size(y)); y];
    Poavx   = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
    power_out = trapz(y*1e-6, Poavx);
    
    % input
    x       = -.5*L;
    xx      = [x*ones(size(y)); y];
    Poavx   = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
    normE   = mphinterp(model, 'ewfd.normE', 'coord', xx);
    power_in= trapz(y*1e-6, Poavx);
    
    [par,~] = fminsearch(@(p) fit_gauss(p, y, normE), ...
        [max(normE) sum(normE.*y)/sum(normE) 1], ooo);
    w_0s(iter)      = par(3);
    n_effs(iter)    = mphglobal(model, 'real(ewfd.neff_1)');
    power_out_sig(iter) = power_out/power_in;
    
    figure(iter); clf;
    plot(y, normE, '-', y, gauss(y, par(1), par(2), par(3)), '--');
    xlabel('x [{\mu}m]');
    ylabel('y [V / m]');
    legend('numerical', 'fitting' );
    error = sum((gauss(y,par(1),par(2),par(3))-normE).^2) ...
        /sum((normE - mean(normE)).^2);
    fprintf('Relative squared error = %2.4e \n', error);

    % collecting dataset
    normEs(iter,:)  = normE;
    fittings(iter,:)= gauss(y, par(1), par(2), par(3));
    errors(iter)    = error;
    
%     if iter == 1 || iter == num_of_ite - 2
%         fprintf('collecting RSE for sigma = %2.2f [um] \n', sigmas(iter));
%         figure(iter); clf;
%         plot(y, normE, '-', y, gauss(y, par(1), par(2), par(3)), '--');
%         xlabel('x [{\mu}m]');
%         ylabel('y [V / m]');
%         legend('numerical', 'fitting' );
%         error = sum((gauss(y,par(1),par(2),par(3))-normE).^2) ...
%             /sum((normE - mean(normE)).^2);
%         fprintf('Relative squared error = %2.4e \n', error);
%         err_ite = err_ite + 1;
%         
%         % collecting dataset
%         normEs(err_ite,:)   = normE;
%         fittings(err_ite,:) = gauss(y, par(1), par(2), par(3));
%         errors(err_ite)     = error;
%     end
end

disp('^--------------------Computation completed--------------------^');
%% saving the dataset
filename    =   ['data/tbg_sigma'  model_spec '.mat'];
save(filename, 'y', 'normE0', 'fitting0', 'w_0', 'n_eff', 'error0', ...
    'sigmas', 'n_effs', 'w_0s', 'theta', 'power_out_sig',           ...
    'normEs', 'fittings', 'errors');