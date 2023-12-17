%% Computing for 3D beam shaping check-up with interference pattern
% By interfering the pump and the required field propagated at the grating
% portion, it should give the required grating profile to produce the
% required field that can propagate at the target plane. This calculation
% is done without considering the pump loss but should be sufficient to
% show that the profile that we get from the BEAM TRACING is appropriate.
%
% Looks at two types:
%   1. plane waves
%   2. tilted Gaussian
% from regular pump without loss and gives both horizontal and vertical
% cross-sections at y = 0 and x = 0, respectively.
%
% The field is defined in define_fields.m

% Pump/beam parameters
lam = 780e-9;
x = (-20:.05:20) * 1e-6;
y = (-20:.05:20) * 1e-6;
z = (-3:.05:3) * 1e-6;
n = [1 1 1];

% Horizontal viewpoint
[E1_pum, E1_tar] = define_fields(x, 0, y, lam, 'gaussian', ...
    n, 'horizontal');
%%
% figure(11); clf;
% pcolor(x, y, abs(squeeze(E1_pum)).^2);
% xlabel('x'); ylabel('y');
% title('|E1_pump|^2, y=0');
% shading flat;
% axis equal;
% % 
% figure(12); clf;
% pcolor(x, y, abs(squeeze(E1_tar)).^2);  
% xlabel('x'); ylabel('y');
% title('|E1_target|^2, y=0');
% shading flat;
% axis equal;

E1_int = abs(E1_tar+E1_pum).^2-abs(E1_tar).^2-abs(E1_pum).^2;

% figure(13); clf;
% pcolor(x, y, squeeze(E1_int));
% xlabel('x'); ylabel('y');
% title('interference, |E1_t+E1_p|^2-|E1_t|^2-|E1_p|^2');
% shading flat
% axis equal

% Vertical viewpoint
[E2_pum, E2_tar] = define_fields(x, z, 0, lam, 'gaussian', ...
    n, 'vertical');

% figure(21); clf;
% pcolor(x, z, abs(squeeze(E2_pum)).^2);
% xlabel('x'); ylabel('z');
% title('|E1_pump|^2, z=0');
% shading flat;
% axis equal;
% 
% figure(22); clf;
% pcolor(x, z, abs(squeeze(E2_tar)).^2);  
% xlabel('x'), ylabel('z')
% title('|E1_target|^2, z=0');
% shading flat;
% axis equal;

E2_int = abs(E2_tar+E2_pum).^2-abs(E2_tar).^2-abs(E2_pum).^2;

% figure(23); clf;
% pcolor(x, z, squeeze(E2_int));
% xlabel('x'); ylabel('z');
% title('interference, |E1_t+E1_p|^2-|E1_t|^2-|E1_p|^2');
% shading flat
% axis equal

figure(30); clf;
ax1 = axes; pcolor(1e6*x, 1e6*y, squeeze(E1_int));
shading interp;
set(gca, 'FontSize', 16);

ax2 = axes; pcolor(1e6*x, 1e6*z, squeeze(E2_int));
shading interp;
ylim([-3 3]);
set(gca, 'FontSize', 16);

ax1.XTick = []; ax1_height = .68; ax2_height = 6*ax1_height/40;
set(ax1, 'Position',[.21 .26 .50 ax1_height]);
set(ax2, 'Position',[.21 .13 .50 ax2_height]);

cb1 = colorbar(ax1,'Position',[.75 .26 .04 ax1_height]);
colorbar(ax2,'Position',[.75 .13 .04 ax2_height]);

xlabel(ax2, 'x/ {\mu}m', 'fontsize', 16);
ylabel(ax1, 'y/ {\mu}m', 'fontsize', 16);
ylabel(cb1, '{\Delta}n_g', 'fontsize', 16);

%% cross section
figure(50); clf;
plot(1e6*x, squeeze(E1_int(401, :)))

figure(51); clf;
plot(1e6*x, squeeze(E2_int(61, :)))