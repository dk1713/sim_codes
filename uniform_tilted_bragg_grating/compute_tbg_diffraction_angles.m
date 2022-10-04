%% Load the grating file
model = mphload('tbg_uniform.mph');
%           
%    -------------power_up--------------
%    |                                 |
% power_in                          power_out
%    |                                 |
%    ------------power_down-------------
%
% mode index: 1.460, 1.460, 1.4599

%% Input Mode
disp('v---------------------Computation started---------------------v')
disp('Computing for width of the input mode')
model.param('par3').set('dn_g',   num2str(0) );
model.study('std1').run

L       = mphglobal(model, 'L');
H       = mphglobal(model, 'H');
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

% Checking for input mode
figure(1)
plot(y, normE, 'b', y, gauss(y, par(1), par(2), par(3)), '--r')
xlabel('x [\mu m]')
ylabel('y [V / m]')
legend('numerical', 'fitting' )
disp('^--------------------Computation completed--------------------^')

%% Filename:
L_g     = mphglobal(model, 'L_g');
H_g     = mphglobal(model, 'H_g');
n_cl    = mphglobal(model, 'n_cl');
n_co    = mphglobal(model, 'n_co');
sigma   = mphglobal(model, 'sigma');
dn      = mphglobal(model, 'dn');

dn_g    = 4e-3;

model_spec = [                      ...
    '_L_',      num2str(L),         ...
    '_H_',      num2str(H),         ...
    '_L_g_',    num2str(L_g),       ...
    '_H_g_',    num2str(H_g),       ...
    '_n_cl_',   num2str(n_cl),      ...
    '_n_co_',   num2str(n_co),      ...
    '_sigma_',  num2str(sigma),     ...
    '_dn_',     num2str(dn),        ...
    '_dng_',    num2str(dn_g),      ...
    ];

%% Init for the parametric sweep
% Setting up the grating length
total_size  = 2^10;
y       = linspace(-.5*H, .5*H, total_size);
x       = .5*L;
xx      = [x*ones(size(y)); y];
model.param('par3').set('dn_g', num2str(dn_g));

%% Transmission against diffraction angles 
disp('v---------------------Computation started---------------------v')
fprintf('Computing for index Case: L_g = %4.3f and varying phi\n', L_g);

% init
phis_deg = linspace(20, 160, 11);
num_of_iter = length(phis_deg);
trans   = zeros(size(phis_deg));
thetas  = zeros(size(phis_deg));

for iter = 1:num_of_iter
    fprintf('Computing for phi = %2.2f [deg]', phis_deg(iter));
    fprintf('    (iteration #%i of %i iterations)\n', iter, num_of_iter);
    
    % Selecting diffraction angle [deg]
    model.param('par2').set('phi', [num2str(phis_deg(iter)), ' [deg]']);
    model.study('std1').run;
    
    % Measuring time-average powerflow in y
    Poavx       = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
    power_out   = trapz(y*1e-6, Poavx);
    
    thetas(iter)  = mphglobal(model, 'theta');
    trans(iter)   = power_out/power_in;  
end
disp('^--------------------Computation completed--------------------^')

%% Save the data generated
filename    =   ['data/tbg_tran_phi'  model_spec '.mat'];
save(filename, 'n_eff', 'w_0', 'power_in', 'thetas', 'dn_g', ...
    'phis_deg', 'trans');