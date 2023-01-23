% 3d interference patterns (holograms)

% interference of a pump beam (Gaussian shape in z, plane wave in x,y) and
% a target beam (plane wave or proper 3D Gaussian)

lam = 780e-9;
k0 = 2*pi/lam;

yzv = [0 0.2 0.4]*1e-6;   % for loop of y or z (cross sections at different positions)
nyzv = length(yzv);

for iyzv=1:nyzv

if (1==0)   % vertical cross section at y=const
    plotsel = 1;
    xv = (-20:0.05:20)*1e-6;
    yv = yzv(iyzv);     %0;     %(-4:0.05:4)*1e-6;
    zv = (-3:0.05:3)*1e-6;
else        % horizontal cross section at z=const
    plotsel = 2;
    xv = (-20:0.05:20)*1e-6;
    yv = (-20:0.05:20)*1e-6;
    zv = yzv(iyzv);     %0.e-6;     %(-3:0.05:3)*1e-6;    
end

[xx,yy,zz] = meshgrid(xv,yv,zv);

% pump field

k1x = k0;
w1 = 2e-6;
w1y = 100*5.8e-6;

E1 = exp(1i*k1x*xx).*exp(-zz.^2/w1^2).*exp(-yy.^2/w1y^2);

figure(101+(iyzv-1)*10)
if plotsel==1
    pcolor(xv,zv,abs(squeeze(E1)').^2)
    xlabel('x'), ylabel('z')
    title(['pump, |E1|^2, y=' num2str(yv)])
else
    pcolor(xv,yv,abs(squeeze(E1)).^2)    
    xlabel('x'), ylabel('y')
    title(['pump, |E1|^2, z=' num2str(zv)])
end
shading flat
axis equal

% target field

if (1==0)   % plane wave
    %al2 = 120 *pi/180;   % propagation angle, 90deg=vertical
    n = [1 1 10];   % direction of propagation
    n = n/norm(n);

    k2x = k0*n(1);
    k2y = k0*n(2);
    k2z = k0*n(3);

    E20 = 0.1;

    E2 = E20*exp(1i*(k2x*xx+k2y*yy+k2z*zz));
end

if (0==0)   % tilted Gaussian
    % al2 = 60 *pi/180;   % propagation angle, 90deg=vertical
    n = [0 0 1];   % direction of propagation
    n = n/norm(n);

    w2 = 2.5e-6;
    zfoc = 50e-6;    % distance to focal point
    E20 = 0.1;

    
    zzr = xx*n(1) + yy*n(2) + zz*n(3);  % z in rotated frame = projection onto prop. direction
    xxr = sqrt(xx.^2+yy.^2+zz.^2 - zzr.^2);     % x=r in rotated frame = the orth. component to z
    zzr = zzr - zfoc;
    
    zR = pi*w2^2/lam;
    w2z = w2*sqrt(1+(zzr/zR).^2);
    eta = atan(zzr/zR);
    %Rz = zzr.*(1+(zR./zzr).^2);
    Rzi = zzr./(zzr.^2 + zR^2);        % inverse of curvature
    
    %E2 = E20 * w2./w2z .* exp(-(xxr./w2z).^2) .* exp(1i*(k0*zzr+k0*xxr./(2*Rz)-eta));
    E2 = E20 * w2./w2z .* exp(-(xxr./w2z).^2) .* exp(1i*(k0*zzr+k0*xxr.^2.*Rzi/2-eta));
end
    
figure(102+(iyzv-1)*10)
if plotsel==1
    pcolor(xv,zv,abs(squeeze(E2)').^2)
    xlabel('x'), ylabel('z')
    title(['target, |E2|^2, y=' num2str(yv)])
else
    pcolor(xv,yv,abs(squeeze(E2)).^2)    
    xlabel('x'), ylabel('y')
    title(['target, |E2|^2, z=' num2str(zv)])
end
shading flat
axis equal

% interference pattern

Et = E1+E2;

figure(103+(iyzv-1)*10)
if plotsel==1
    pcolor(xv,zv,abs(squeeze(Et)').^2)
    xlabel('x'), ylabel('z')
    title(['interference, |E1+E2|^2, y=' num2str(yv)])
else
    pcolor(xv,yv,abs(squeeze(Et)).^2)    
    % pcolor(xv,yv,(abs(squeeze(Et)).^2-abs(squeeze(E1)).^2-abs(squeeze(E2)).^2)./abs(squeeze(E1)).^2)    
    %   note on above:  this is meant to be ~dng, but shouldn't be
    %                   normalised to pump field at z, only to the pump at centre of
    %                   waveguide, or power integrated over z
    xlabel('x'), ylabel('y')
    title(['interference, |E1+E2|^2, z=' num2str(zv)])
end
shading flat
axis equal

% interference term only

Et2 = abs(E1+E2).^2-abs(E1).^2-abs(E2).^2;

figure(104+(iyzv-1)*10)
if plotsel==1
    pcolor(xv,zv,squeeze(Et2)')
    xlabel('x'), ylabel('z')
    title(['interference, |E1+E2|^2-|E1|^2-|E2|^2, y=' num2str(yv)])
else
    pcolor(xv,yv,squeeze(Et2))    
    % pcolor(xv,yv,(abs(squeeze(Et)).^2-abs(squeeze(E1)).^2-abs(squeeze(E2)).^2)./abs(squeeze(E1)).^2)    
    %   note on above:  this is meant to be ~dng, but shouldn't be
    %                   normalised to pump field at z, only to the pump at centre of
    %                   waveguide, or power integrated over z
    xlabel('x'), ylabel('y')
    title(['interference, |E1+E2|^2-|E1|^2-|E2|^2, z=' num2str(zv)])
end
shading flat
axis equal

end
