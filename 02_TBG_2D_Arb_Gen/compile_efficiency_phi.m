load('efficiency_gauss.mat')
x2 = x3;
s2 = s3;
load('efficiency_super.mat')
%%
figure(1)
plot(x2 + 90, s2, 'r', x3 + 90, s3, 'b')
xlabel('scattering angle, \phi / [deg]')
ylabel('reflectance efficiencies [%]')
legend('Gaussian', 'super-Gaussian')