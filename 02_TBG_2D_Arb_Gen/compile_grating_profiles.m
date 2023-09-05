% load('realistic_super.mat')
load('realistic_gauss.mat')
%% parameters
lam     = 780e-9;
phi     = 30*pi/180;
n_eff   = 1.4635;

%%
figure(10); clf;
colororder(["#8040E6";"#1AA640"])
yyaxis left
plot(1e3*x, dn_gs, 'linewidth', 1);
xlabel('x/ mm');
ylabel('index modulation, {\Delta}n_g')

yyaxis right
plot(1e3*x, 1e9*period, 'linewidth', 1 );
yline(1e9*lam/(1*cos(phi+.5*pi) + n_eff), '--k', 'constant k');
ylabel('grating period, {\Lambda}(x) / [nm]')

xlim([-3, 3]);