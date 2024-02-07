%% plotting mode shape against the Gaussian fitting
% This is a simple plotting script. It takes in fitting_[sigma].mat
% sigma values are: 0.5um, 2um and 5um

sigma = [0.5, 2, 5];

for ite = 1:length(sigma)
    load(['data/fitting' num2str(sigma(ite)) '.mat']);
    figure(ite); clf;
    plot(y, normE, y, fitting, '--', 'LineWidth', 2);
    xlabel('y / [{\mu}m]');
    ylabel('electric field norm / [Vm^{-1}]');
    legend('numerical', 'fitting', 'Location', 'best');
    set(gca, 'FontSize', 16);
end