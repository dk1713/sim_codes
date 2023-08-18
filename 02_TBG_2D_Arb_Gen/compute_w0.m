%% Compute the width of the guide mode
% This is needed to set the w0 inside the COMSOL. May need to save the
% COMSOL file manually since saving in the script sometimes doesn't work.

% Run Comsol
model = mphload('grating_design_arbitrary_v5.3.mph');
model.param('par3').set('control', num2str(0));
model.study('std1').run

L       = mphglobal(model, 'L');
H       = mphglobal(model, 'H');

y       = linspace(-.5*H, .5*H, 500);
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
model.param('par3').set( 'w_0', [num2str(par(3)), '[um]'] );
model.param('par3').set('control', num2str(1));
% model.save('grating_design_arbitrary_v5.3.mph');
%%
figure(1)
plot(y, normE, 'b', y, gauss(y, par(1), par(2), par(3)), '--r')
xlabel('x [\mu m]')
ylabel('y [V / m]')
legend('numerical', 'fitting' )