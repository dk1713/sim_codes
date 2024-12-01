% Fresnel coefficients
n1 = 1;
n2 = 1.44;

theta = linspace(0,60,500)*pi/180;

R_s = 1./(2*cos(theta).^2).^2;
R_p = R_s .* cos(2*theta).^2;

figure(1)
plot(180*theta/pi, R_s, 180*theta/pi, R_p)
xlabel('tilt angle, \theta / [deg]');
ylabel('Scattering intensity');
legend('s-pol','p-pol','location','northwest')