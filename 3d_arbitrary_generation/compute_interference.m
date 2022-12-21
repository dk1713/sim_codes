%% Computeing for grating profile
% Takes in the data from compute_3d_arb_beam.m to generate the grating
% profile needed to generate arbitrary beam specified in that script.

% Indices [1]
n_air   = 1;
n_clad  = 1.4555;
n_core  = 1.4608;
n_eff   = 1.4635;
dn_g    = 3e-3;

% Heights of the layers [m]
h_clad  = 15e-6;
h_core  = 5e-6;
sigma   = .5*h_core;
lambda  = 780e-9;

load('3d_gauss.mat');
beta = 2*pi*n_eff/lambda;

% init memo allo
z_g     = linspace(-sigma, sigma, 2^10);
n_grat  = zeros([size(EE_grat), length(z_g)]);

% input mode
E_in      = exp(1i*beta*xx_g).*power_ratio';
%% Checking the plot
figure(1)
pcolor(xx_g, yy_g, abs(E_in'))
shading flat
xlabel('x'), ylabel('y')

figure(2)
pcolor(xx_g, yy_g, abs(EE_grat).^2)g
shading flat
xlabel('x'), ylabel('y')

figure(3)
pcolor(xx_g, yy_g, abs(E_in' + EE_grat).^2 - abs(E_in').^2- abs(EE_grat).^2)
shading flat
xlabel('x'), ylabel('y')
%%
% for i = 1:length(z_g)
%     n_grat(:,:,i) = abs(E_in'*exp(-(z_g(i)/sigma)^2) );
% end
% 
% %%
% figure(10)
% for i = 1:10:length(z_g)
%     filename = 'animation.gif';
%     [~,h] = contourf(1e6*xx_g, 1e6*yy_g, n_grat(:,:,i));
%     h.EdgeColor = 'none';
%     colorbar;
%     axis equal
%     % Updating the title
%     title(sprintf('Refractive index\n z = %2.4f [um]', 1e6*z_g(i)),...
%     'Interpreter','Latex');
%     
%     % Delay
%     pause(0.1)
%     % Saving the figure
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
%         'DelayTime',1);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append',...
%         'DelayTime',1);
%     end
% end

%%
% figure(100)
% hold on
% for i = 1:2:length(z_g)
% %     surf(1e3*xx_g, 1e3*yy_g, n_grat(:,:,i), 'EdgeColor', 'none');
% %     contour(1e3*xx_g, 1e3*yy_g, n_grat(:,:,i));
%     [~,h] = contourf(1e6*xx_g, 1e6*yy_g, n_grat(:,:,i));
%     h.EdgeColor = 'none';
%     h.ContourZLevel = 1e6*z_g(i);
%     axis equal
% end
% hold off
% colorbar;
% view(3);
% 
% %%
% [xg, zg] = meshgrid(xx_g(1,:), z_g);
% temp = zeros(size(xg));
% 
% for i = 1:size(xg,1)
%     for j = 1:size(xg,2)
%         temp(i,j) = n_grat(128,j,i);
%     end
% end
% 
% figure(202)
% contourf(1e6*xg, 1e6*zg, temp, 'EdgeColor','none');
% axis equal
% colorbar;