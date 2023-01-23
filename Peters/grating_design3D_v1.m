% design of a tilted grating for a 3D arbitrarily shaped output beam

% assumptions:
% - target field scalar, from scattering formula match |E| in plane (without the vertical component) to target
% - grating tilt angle is the same everywhere (not correct? local change of
%   grating leads to change of effective tilt angle?)
% - neglect diffraction of pump beam

global lam c k0 neff beta w0 sigma dng

% parameters

lam = 780e-9;           % wavelength
k0 = 2*pi/lam;

c = 3e8;

neff = 1;  %NOT IMPLEMENTED EVERYWHERE           % effective index of fundamental mode
beta = neff*k0;         % propagation constant of mode
w0 = 2e-6;              % waist of fundamental mode (vertical)
sigma = 2e-6;           % waist of refractive index profile (vertical)
dng = 1;                % grating index contrast (=1 for design calculation; realistic <=5e-3 ?)


% define target beam (scalar), propagate to grating plane

    % grid in grating plane: horizontal cross section at z=0
    xv = (-20:0.05:20)*1e-6;
    yv = (-20:0.05:20)*1e-6;

    [xx,yy] = meshgrid(xv,yv);

    % target field

    if (1==0)   % plane wave
        %al2 = 120 *pi/180;   % propagation angle, 90deg=vertical
        ntar = [1 1 1];   % direction of propagation
        ntar = ntar/norm(ntar);

        k2x = k0*ntar(1);
        k2y = k0*ntar(2);
        k2z = k0*ntar(3);

        E20 = 0.1;      % target field amplitude - THIS SHOULD NOT HAVE ANY EFFECT

        E2 = E20*exp(1i*(k2x*xx+k2y*yy));
    end

    if (0==0)   % tilted Gaussian
        % al2 = 60 *pi/180;   % propagation angle, 90deg=vertical
        ntar = [1 1 1];   % direction of propagation
        ntar = ntar/norm(ntar);

        w2 = 2.5e-6;
        zfoc = 50e-6;    % distance to focal point
        E20 = 0.1;       % target field amplitude - THIS SHOULD NOT HAVE ANY EFFECT

        zzr = xx*ntar(1) + yy*ntar(2);  % z in rotated frame = projection onto prop. direction
        xxr = sqrt(xx.^2+yy.^2 - zzr.^2);     % x=r in rotated frame = the orth. component to z
        zzr = zzr - zfoc;

        zR = pi*w2^2/lam;
        w2z = w2*sqrt(1+(zzr/zR).^2);
        eta = atan(zzr/zR);
        %Rz = zzr.*(1+(zR./zzr).^2);
        Rzi = zzr./(zzr.^2 + zR^2);        % inverse of curvature

        %E2 = E20 * w2./w2z .* exp(-(xxr./w2z).^2) .* exp(1i*(k0*zzr+k0*xxr./(2*Rz)-eta));
        E2 = E20 * w2./w2z .* exp(-(xxr./w2z).^2) .* exp(1i*(k0*zzr+k0*xxr.^2.*Rzi/2-eta));
    end
    
    figure(1)
    pcolor(xv,yv,abs(E2).^2)    
    xlabel('x'), ylabel('y')
    title(['target, |E2|^2, z=0'])
    shading flat
    axis equal
    
    figure(2)
    pcolor(xv,yv,angle(E2))    
    xlabel('x'), ylabel('y')
    title(['target, phase(E2), z=0'])
    shading flat
    axis equal

% calculate central grating properties (tilt, rotation, period)

[lamgrat0,alphagrat0,alphatilt0] = grating_angles_3D_f2(ntar(1),ntar(2),ntar(3));
lamgrat0 = lamgrat0 * lam;

% at each position, calculate local target propagation direction, local grating period

    % local propagation direction of target field, extracted from E-grid
    kx = -1i * (E2(:,3:end)-E2(:,1:end-2)) / (2*(xv(2)-xv(1)));
    kx = kx ./ E2(:,2:end-1);
    kx = [kx(:,1) kx kx(:,end)];
    kx = real(kx);

    ky = -1i * (E2(3:end,:)-E2(1:end-2,:)) / (2*(yv(2)-yv(1)));
    ky = ky ./ E2(2:end-1,:);
    ky = [ky(1,:); ky; ky(end,:)];
    ky = real(ky);
    
    kz = sqrt(k0^2-kx.^2-ky.^2);

    [lamgrat,alphagrat,alphatilt] = grating_angles_3D_f2(kx,ky,kz);
    lamgrat = lamgrat * lam;

    figure(3)
    subplot(221)
    pcolor(xv,yv,kx)    
    xlabel('x'), ylabel('y')
    title('kx of local target field')
    shading flat
    axis equal
    colorbar
    
    subplot(222)
    pcolor(xv,yv,ky)    
    xlabel('x'), ylabel('y')
    title('ky of local target field')
    shading flat
    axis equal
    colorbar

    subplot(223)
    pcolor(xv,yv,kz)    
    xlabel('x'), ylabel('y')
    title('kz of local target field')
    shading flat
    axis equal
    colorbar

    figure(4)
    subplot(221)
    pcolor(xv,yv,lamgrat)    
    xlabel('x'), ylabel('y')
    title('Local grating period, lamgrat')
    shading flat
    axis equal
    colorbar
    
    subplot(222)
    pcolor(xv,yv,alphagrat*180/pi)    
    xlabel('x'), ylabel('y')
    title('Local grating direction, alphagrat (deg)')
    shading flat
    axis equal
    colorbar

    subplot(223)
    pcolor(xv,yv,alphatilt*180/pi)    
    xlabel('x'), ylabel('y')
    title('Local optimum grating tilt, alphatilt (deg)')
    shading flat
    axis equal
    colorbar

% calculate the modified local grating tilt, given the local grating
% direction if it is written with the central tilt angle

    % for initial development, may just assume that the effective tilt
    % angle is the same as the central tilt angle (i.e. omit this step)

    
% calculate scattering rate for E-field and power (apart from dng factor)

pol = 1;            % =1 horizontal input polarisation; =2 vertical input polarisation
ky = 0*lamgrat;     % transverse prop. constant of pump

[al1,Ex1,Ey1,Ez1] = scatteringrate3_f(lamgrat,alphagrat,alphatilt,ky,pol);

E1norm = sqrt(Ex1.^2+Ey1.^2+Ez1.^2);

    figure(5)
    pcolor(xv,yv,al1)    
    xlabel('x'), ylabel('y')
    title('Scattering efficiency al (1/m) for dng=1')
    shading flat
    axis equal
    colorbar

    figure(6)
    subplot(221)
    pcolor(xv,yv,Ex1./E1norm)    
    xlabel('x'), ylabel('y')
    title('Scattered field pol.: Ex (norm.)')
    shading flat
    axis equal
    colorbar

    subplot(222)
    pcolor(xv,yv,Ey1./E1norm)    
    xlabel('x'), ylabel('y')
    title('Scattered field pol.: Ey (norm.)')
    shading flat
    axis equal
    colorbar

    subplot(223)
    pcolor(xv,yv,Ez1./E1norm)    
    xlabel('x'), ylabel('y')
    title('Scattered field pol.: Ez (norm.)')
    shading flat
    axis equal
    colorbar

% from this, extract dng (to get correct E) - for no pump depletion, this is all
%   NOTE: there would be different possibilities, because the target field is
%   scalar, but the scattered field has polarisation, e.g. trying to match
%   modulus of E-field, match one polarisation, ...
%   NOTE: can also include pump profile here (e.g. Gaussian in transverse
%   direction)

dng1 = abs(E2)./E1norm;
dng1 = dng1/(max(dng1(:)));

    figure(7)
    subplot(221)
    pcolor(xv,yv,dng1)    
    xlabel('x'), ylabel('y')
    title('Grating dng (norm.), no pump depletion')
    shading flat
    axis equal
    colorbar

    subplot(222)
    contour(xv,yv,dng1)    
    xlabel('x'), ylabel('y')
    title('Grating dng (norm.), no pump depletion')
    shading flat
    axis equal
    colorbar

% to include pump depletion, integrate the grating "loss" and scale dng
% accordingly larger with grating propagation distance to maintain correct
% E field

loss = cumsum(dng1.^2 .* al1,2) * (xv(2)-xv(1));

selpump = 3;
if (selpump==1)    % case 1: flat pump at input P=1, choose max dng for pump depletion
    dngbar_max = sqrt(1/max(loss(:,end)));

    dngbar = dngbar_max*0.99;
    P = 1 + 0*xx;
    P = P - dngbar^2*loss;
    dng_full = real(dng1.*dngbar./sqrt(P));
end
if (selpump==2)    % case 2: flat pump at input P=1, choose max dng for fractionial pump depletion
    dngbar_max = sqrt(1/max(loss(:,end)));

    dngbar = dngbar_max*sqrt(0.01);
    P = 1 + 0*xx;
    P = P - dngbar^2*loss;
    dng_full = real(dng1.*dngbar./sqrt(P));
end
if (selpump==3)    % case 3: optimised pump at input, choose max dng for 1/2 pump depletion
    dngbar_max = sqrt(1/max(loss(:,end)));

    dngbar = dngbar_max*0.99;
    P = (dngbar_max^2*loss(:,end))*ones(1,length(xv));
    P = P - dngbar^2*loss;
    dng_full = real(dng1.*dngbar./sqrt(P));
end

    figure(8)
    subplot(221)
    pcolor(xv,yv,P)    
    xlabel('x'), ylabel('y')
    title('Pump power')
    shading flat
    axis equal
    colorbar

    subplot(223)
    pcolor(xv,yv,dng_full)   
    xlabel('x'), ylabel('y')
    title('Grating dng, pump depletion')
    shading flat
    axis equal
    colorbar

    subplot(224)
    contour(xv,yv,dng_full)   
    xlabel('x'), ylabel('y')
    title('Grating dng, pump depletion')
    shading flat
    axis equal
    colorbar

% combine this grating strength with phase information to produce schematic
% of the full grating (to compare with the holographic results)

xxplot = xx;
yyplot = yy;
zzplot = 0*xx;

%alphagrat0
phi = angle(E2.*exp(-1i*k0*(ntar(1)*xxplot+ntar(2)*yyplot+ntar(3)*zzplot)));

ng = dng_full.*sin(2*pi/lamgrat0*(xxplot*cos(alphagrat0)+yyplot*sin(alphagrat0)-zzplot*tan(alphatilt0)) - phi);

    figure(9)
    %subplot(221)
    pcolor(xv,yv,ng)    
    xlabel('x'), ylabel('y')
    title('Full grating ng')
    shading flat
    axis equal
    colorbar

xvp = xx(1,:);
zvp = (-3:0.05:3)*1e-6;

[xxplot,zzplot] = meshgrid(xvp,zvp);
yyplot = 0*xxplot;

phip = phi((length(xx(:,1))+1)/2,:);     % phase of target field at y=0
phip = ones(length(zvp),1) * phip;

dng_fullp = dng_full((length(xx(:,1))+1)/2,:);     % dng at y=0
dng_fullp = ones(length(zvp),1) * dng_fullp;

ngp = dng_fullp.*sin(2*pi/lamgrat0*(xxplot*cos(alphagrat0)+yyplot*sin(alphagrat0)-zzplot*tan(alphatilt0)) - phip).*exp(-(zzplot/sigma).^2);

    figure(10)
    %subplot(221)
    pcolor(xvp,zvp,ngp)    
    xlabel('x'), ylabel('z')
    title('Full grating ng')
    shading flat
    axis equal
    colorbar

% finally, calculate the scattered field including phase information,
% propagate to target position, check result



