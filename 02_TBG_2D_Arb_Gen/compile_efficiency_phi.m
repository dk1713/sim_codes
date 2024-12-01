%%
load('data/effi_dist=8mm/efficiency_super.mat')
x2 = x;
s2 = s;
load('data/effi_dist=8mm/efficiency_gauss.mat')
%%
% load('data/effi_dist=5mm/efficiency_super.mat')
% x2 = x;
% s2 = s;
% load('data/effi_dist=5mm/efficiency_gauss.mat')
%%
figure();clf;
plot(x + 90, s, 'r', x2 + 90, s2, 'b', 'linewidth', 2);
xlabel('scattering angle, \phi / [deg]')
ylabel('diffraction efficiencies [%]')
legend('Gaussian', 'super-Gaussian')