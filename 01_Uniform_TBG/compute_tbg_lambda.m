%% Load the grating file
model = mphload('tbg_constant_dng.mph');
%           
%    -------------power_up--------------
%    |                                 |
% power_in                          power_out
%    |                                 |
%    ------------power_down-------------
%

%% Taking up parameter data from .mph
% for the filename labelling
L       = mphglobal(model, 'L');
H       = mphglobal(model, 'H');
L_g     = mphglobal(model, 'L_g');
H_g     = mphglobal(model, 'H_g');
n_cl    = mphglobal(model, 'n_cl');
n_co    = mphglobal(model, 'n_co');
sigma   = mphglobal(model, 'sigma');
dn      = mphglobal(model, 'dn');

% Introduce lower dng
dng = 3e-3;
model.param('par3').set( 'dn_g',   num2str(dng) );

% Selecting scattering angle [deg]
phi = [60, 90, 120];
for iter = 1:length(phi)
    phi_deg     = phi(iter);
    model.param('par2').set( 'phi',     [num2str(phi_deg),      ' [deg]'] );
    
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
        '_dng_',    num2str(dng),   ...
        ];

    %% Init for three parametric sweep
    total_size  = 2^9;

    % function handles
    gauss       =   @(x, a, b, c)...
        a*exp(-((x-b)/c).^2);
    fit_gauss   =   @(p, x, data)...
        sqrt(mean( (gauss(x,p(1),p(2),p(3)) - data).^2  ));

    %% lambdas efficiency studies
    disp('v---------------------Computation started---------------------v');
    fprintf('Computing for index Case: dng = %2.3e and varying lambda\n', dng);

    lambdas     = linspace(740, 820, 9);
    num_of_iterations = length(lambdas);
    power_loss_lambda = ones(size(lambdas));
    n_effs      = zeros(size(lambdas));
    w_0s        = zeros(size(lambdas));
    phis        = zeros(size(lambdas));

    % Measurement line
    y = linspace(-.5*H, .5*H, total_size);
    
    for iter2 = 1:num_of_iterations
        fprintf('Computing for lambda = %4.1f [nm]', lambdas(iter2));
        fprintf('    (iteration #%i of %i iterations)\n', ...
            iter2, num_of_iterations);

        model.param('par2').set('lambda0',  [num2str(lambdas(iter2)), ' [nm]']);
        model.study('std1').run;
        
        % power in
        xx      = [-.5*L*ones(size(y)); y];  
        normE   = mphinterp(model, 'ewfd.normE', 'coord', xx);
        Poavx   = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
        power_in = trapz(y*1e-6, Poavx);

        % mode width
        ooo = optimset('TolX',1e-9);
        [par,~] = fminsearch(@(p) fit_gauss(p, y, normE), ...
                [max(normE) sum(normE.*y)/sum(normE) 1], ooo);
        w_0s(iter2) = par(3);

        % n_eff
        n_effs(iter2) = mphglobal(model, 'real(ewfd.neff_1)');

        xx      = [.5*L*ones(size(y)); y];       
        Poavx = mphinterp(model, 'ewfd.Poavx', 'coord', xx);
        power_loss_lambda(iter2) = trapz(y*1e-6, Poavx)/power_in;
        
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
        
        phis(iter2)= pi/2 + atan(Sx/Sy);
    end

    disp('^--------------------Computation completed--------------------^');
    %%
    filename    =   ['data/tbg_lambda'  model_spec '.mat'];
    save(filename, 'n_effs', 'w_0s', 'lambdas', 'phis', 'power_loss_lambda');
end