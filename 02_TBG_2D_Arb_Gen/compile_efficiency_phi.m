load('efficiency_super.mat')
x2 = x;
s2 = s;
load('efficiency_gauss.mat')
%%
figure()
plot(x + 90, s, 'r', x2 + 90, s2, 'b')
xlabel('scattering angle, \phi / [deg]')
ylabel('reflectance efficiencies [%]')
legend('Gaussian', 'super-Gaussian')