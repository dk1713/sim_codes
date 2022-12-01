% grating periods, directions, tilt angles for 3D beam shaping

% x = pump forward direction; y = transverse (in-plane), z = vertical (out of plane)

lam = 1;            % lengths normalised to wavelength
k0 = 2*pi/lam;      % k in material, assumed equal to beta

xtar = 0; -1:0.01:1;   %0;
ytar = 1;           %-1:0.01:1;
ztar = 1;

koutx = k0*xtar./sqrt(xtar.^2+ytar.^2+ztar.^2);
kouty = k0*ytar./sqrt(xtar.^2+ytar.^2+ztar.^2);
koutz = k0*ztar./sqrt(xtar.^2+ytar.^2+ztar.^2);

kinx = k0;
kiny = 0;
kinz = 0;

kgratx = kinx-koutx;    % grating K
kgraty = kiny-kouty;

kgratn = sqrt(kgratx.^2+kgraty.^2); % norm of grating K
lamgrat = 2*pi./kgratn;             % grating period
alphagrat = atan(kgraty./kgratx);   % grating direction (angle)

% normal to the grating planes

ngratx = koutx-kinx;
ngraty = kouty-kiny;
ngratz = koutz-kinz;
nn = sqrt(ngratx.^2+ngraty.^2+ngratz.^2);
ngratx = ngratx./nn;
ngraty = ngraty./nn;
ngratz = ngratz./nn;
alphatilt = pi/2-acos(ngratz);   % grating tilt angle w.r.t. vertical (=0 vertical, =90 deg horizontal planes)


figure(1)
subplot(221)
plot(xtar./ztar,lamgrat)
xlabel('xtar/ztar')
ylabel('lam grating/lam')
grid

subplot(222)
plot(xtar./ztar,kgratn/k0)
xlabel('xtar/ztar')
ylabel('K grating/beta')
grid

subplot(223)
plot(xtar./ztar,alphagrat*180/pi)
xlabel('xtar/ztar')
ylabel('grating direction (deg)')
grid

subplot(224)
plot(xtar./ztar,alphatilt*180/pi)
xlabel('xtar/ztar')
ylabel('tilt angle to vertical (deg)')
grid

lamgrat

alphatilt*180/pi

alphagrat*180/pi
