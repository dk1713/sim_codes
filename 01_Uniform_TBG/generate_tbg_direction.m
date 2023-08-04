%% Generation of the directional output dataset from COMSOL simulations
% Preset COMSOL file: 'tbg_constant_dng.mph'
model = mphload('tbg_constant_dng.mph');
% This script generates the numerical dataset from the COMSOL simulation.
% Note that this takes time and might need to check if preset COMSOL file
% works beforehand. (see generate_tbg.m for details on geometry)

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

%% Input Mode
disp('v---------------------Computation started---------------------v');
disp('Computing for width of the input mode');
% setting no dng s.t. there is no back reflection. 
model.param('par3').set( 'dn_g',   num2str(0) );
model.study('std1').run

% n_eff
n_eff = mphglobal(model, 'real(ewfd.neff_1)');

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

ooo = optimset('TolX', 1e-9);
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

%% Init
total_size      = 2^9;
lambdas         = linspace(740, 820, 9);
num_of_ite      = length(lambdas);
power_in_lam    = zeros(size(lambdas));
power_out_lam   = zeros(size(lambdas));
n_effs      	= zeros(size(lambdas));
w_0s            = zeros(size(lambdas));
phis            = zeros(size(lambdas));

% setting higher dng
dng = 3e-3;
model.param('par3').set( 'dn_g',   num2str(dng) );

% Selecting scattering angle [deg]
phi = [60, 90, 120];
for iter1 = 1:length(phi)
    % Selecting scattering angle [deg]
    phi_deg     = phi(iter1);
    model.param('par2').set( 'phi',     [num2str(phi_deg),     ' [deg]'] );

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

    % efficiency calculation
    disp('v---------------------Computation started---------------------v');
    fprintf('Computing for index Case: dn = %2.3e and varying lambda\n', dn);

    for iter2 = 1:num_of_ite
        fprintf('Computing for lambda = %4.1f [nm]', lambdas(iter2));
        fprintf('    (iteration #%i of %i iterations)\n', ...
            iter2, num_of_ite);

        model.param('par2').set('lambda0', [num2str(lambdas(iter2)), ' [nm]']);
        model.study('std1').run;
        
        % n_eff
        n_effs(iter2) = mphglobal(model, 'real(ewfd.neff_1)');
        
        % Measuring enorm
        y       = linspace(-.75*H_g, .75*H_g, total_size);
        x       = -.5*L;
        xx      = [x*ones(size(y)); y];
        normE   = mphinterp(model, 'ewfd.normE', 'coord', xx);
        Poavx   = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
        power_in_lam(iter2)= trapz(y*1e-6, Poavx);
        
        % fitting
        ooo = optimset('TolX',1e-9);
        [par,~] = fminsearch(@(p) fit_gauss(p, y, normE), ...
                [max(normE) sum(normE.*y)/sum(normE) 1], ooo);
        w_0s(iter2) = par(3);

        % Measuring time-average powerflow in y
        x       = .5*L;
        xx      = [x*ones(size(y)); y];
        Poavx   = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
        power_out_lam(iter2) = trapz(y*1e-6, Poavx);
        
        % beam propagation angle
        x       = linspace(-.5*L, .5*L, total_size);
        xx      = [x; (.5*H -2)*ones(size(x))];
        Ez      = mphinterp(model, 'ewfd.Ez', 'coord', xx);
        Ek      = fftshift( fft( fftshift(Ez) ) );
        
        k0      = 2*pi*n_cl/(lambdas(iter2)*1e-9);
        k       = linspace(-pi, pi, length(x) )./(x(1) - x(2))/1e-6;
        ky      = real( sqrt(k0^2 - k.^2) );

        % find angle of beam propagation
        ii      = find(ky>0);
        Sx      = sum(abs(Ek(ii)).^2 .* k(ii));
        Sy      = sum(abs(Ek(ii)).^2 .* ky(ii));
        
        phis(iter2) = pi/2 + atan(Sx/Sy);
    end
    disp('^--------------------Computation completed--------------------^');

    % saving the dataset
    filename    =   ['data/tbg_direction'  model_spec '.mat'];
    save(filename, 'n_eff', 'w_0', 'power_in', 'theta', ...
        'n_effs', 'w_0s', 'lambdas', 'phis', ...
        'power_in_lam', 'power_out_lam');
end