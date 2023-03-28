function [al,Ex,Ey,Ez]=compute_scattering(Lambda,al_grat,theta,w0,sigma,beta,dng,neff,kz,pol)
% formula of power loss from tilted grating
% input: arrays (all of same size) of: 
%       local grating wavelength Lambda (m), 
%       grating direction algrat (rad), 
%       tilt angle (rad), 
%       transverse pump propagation constant kz (1/m)
%   one value for   pol=1 (input horizontal polarisation) or
%                   pol=2 (input vertical polarisation)
% output: corresponding arrays:
%       al = energy loss rate (1/m)
%       Ex, Ey,Ez = scattered electrid field
% NOTE: in output x=pump propagation direction, y=in-plane transverse, z=vertical
%       but internally in this function y=vertical, z=transverse

plotfig = 0;    % for testing: =1 for detailed figures, =0 for no figures

w = 1/sqrt(1/w0^2+1/sigma^2);   % effective width of mode and grating

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% pump at an angle inside slab %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%kz = beta*(-0.3:0.01:0.3);     % transverse k
kx = sqrt(beta^2-kz.^2);        % longitudinal k
al_prop = asin(kz/beta);       % propagation angle


% incident angle relative to grating planes
% ngrat = [-cos(theta), sin(theta), 0];     
% normal to grating plane (x,y,z), then rotate in (x,z) plane by al_grat
% nprop = [cos(al_prop), 0, sin(al_prop)];  % prop. direction

ngratx = -cos(theta).*cos(al_grat);   
ngraty = sin(theta);  
ngratz = -cos(theta).*sin(al_grat);   

npropx = cos(al_prop);  npropy = 0*npropx;  npropz = sin(al_prop);
    
al_inc = acos(-ngratx.*npropx -ngraty.*npropy -npropz.*ngratz);
% figure(112)
% pcolor(180*al_prop/pi)    
% xlabel('x'), ylabel('y')
% shading flat
% axis equal
% colorbar

if plotfig
    figure(12)
    plot(al_prop*180/pi,al_inc*180/pi)
    xlabel('propagation angle')
    ylabel('incident angle rel. grating plane')
end

% define polarisation vectors of s and p, for input (pump) and output (Bragg scattered beam):

% define incident plane, i.e. its normal vector (normal to propagation and normal of grating)
% ninc = ngrat x nprop

nincx = ngraty.*npropz - ngratz.*npropy;
nincy = ngratz.*npropx - ngratx.*npropz;
nincz = ngratx.*npropy - ngraty.*npropx;

nn = sqrt(nincx.^2+nincy.^2+nincz.^2);
nincx = nincx./nn;
nincy = nincy./nn;
nincz = nincz./nn;

% input s polarisation (same as normal to incident plane)
ninsx = nincx; 
ninsy = nincy; 
ninsz = nincz;

% input p pol is normal to propagation direction and s-pol vector, ninp ~ nins x nprop
ninpx = ninsy.*npropz - ninsz.*npropy;
ninpy = ninsz.*npropx - ninsx.*npropz;
ninpz = ninsx.*npropy - ninsy.*npropx;

% propagation direction of scattered beam nscat = nprop - 2*ngrat*(nprop.ngrat)
qprojection = npropx.*ngratx + npropy.*ngraty + npropz.*ngratz;

nscatx = npropx - 2*ngratx.*qprojection;
nscaty = npropy - 2*ngraty.*qprojection;
nscatz = npropz - 2*ngratz.*qprojection;

if plotfig
    figure(30)
    plot(al_prop*180/pi,nscatx,al_prop*180/pi,nscaty,al_prop*180/pi,nscatz)
    xlabel('propagation angle')
    ylabel('x,y,z projections')
    legend('x component','y component','z component')
    title('propagation direction of beam reflected by single grating plane')
end

% output s polarisation: same as input s polarisation

noutsx = nincx; noutsy = nincy; noutsz = nincz;

% p pol is normal to new propagation direction and s-pol vector, noutp ~ nouts x nscat

noutpx = noutsy.*nscatz - noutsz.*nscaty;
noutpy = noutsz.*nscatx - noutsx.*nscatz;
noutpz = noutsx.*nscaty - noutsy.*nscatx;

if plotfig
    figure(31)
    plot(al_prop*180/pi,ninpx,al_prop*180/pi,ninpy,al_prop*180/pi,ninpz)
    xlabel('propagation angle')
    ylabel('x,y,z projections')
    legend('x component','y component','z component')
    title('p-polarised input, x,y,z components')

    figure(32)
    plot(al_prop*180/pi,ninsx,al_prop*180/pi,ninsy,al_prop*180/pi,ninsz)
    xlabel('propagation angle')
    ylabel('x,y,z projections')
    legend('x component','y component','z component')
    title('s-polarised input, x,y,z components')

    figure(33)
    plot(al_prop*180/pi,noutpx,al_prop*180/pi,noutpy,al_prop*180/pi,noutpz)
    xlabel('propagation angle')
    ylabel('x,y,z projections')
    legend('x component','y component','z component')
    title('p-polarised output, x,y,z components')

    figure(34)
    plot(al_prop*180/pi,noutsx,al_prop*180/pi,noutsy,al_prop*180/pi,noutsz)
    xlabel('propagation angle')
    ylabel('x,y,z projections')
    legend('x component','y component','z component')
    title('s-polarised output, x,y,z components')
end

% input vertical (y) polarisation

% direction of polarisation
nex = 0;
ney = 1;
nez = 0;

% project of polarisation on s and p

pys = nex.*ninsx + ney.*ninsy + nez.*ninsz;
pyp = nex.*ninpx + ney.*ninpy + nez.*ninpz;

% input horizontal (x-z) polarisation

%ne = [-sin(al_prop), 0, cos(al_prop)];     % direction of polarisation, orthogonal to nprop
nex = -sin(al_prop);
ney = 0;
nez = cos(al_prop);

% project of polarisation ons and p

pxs = nex.*ninsx + ney.*ninsy + nez.*ninsz;
pxp = nex.*ninpx + ney.*ninpy + nez.*ninpz;

% figure(13)
% plot(al_prop*180/pi,pxs,al_prop*180/pi,pxp,al_prop*180/pi,pys,al_prop*180/pi,pyp)
% xlabel('propagation angle')
% ylabel('power in s and p polarisation')
% legend('z pol, s component','z pol, p component','y pol, s component','y pol, p component')

if plotfig
    figure(13)
    plot(al_prop*180/pi,pxs.^2,al_prop*180/pi,pys.^2)
    xlabel('propagation angle')
    ylabel('power in s polarisation')
    legend('z pol input','y pol input')
end

% output propagation direction and angle (from diffraction)

Kgrat = 2*pi./Lambda;   % grating "propagation" constant

kdiffx = kx - Kgrat.*cos(al_grat);      % diffracted beam
kdiffz = kz - Kgrat.*sin(al_grat);
kdiffy = sqrt(beta^2-kdiffx.^2-kdiffz.^2);

phivar3d = real(asin(kdiffy/beta));

if plotfig
    figure(14)
    plot(al_prop*180/pi,phivar3d*180/pi)
    xlabel('propagation angle')
    ylabel('diffracted beam output angle')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scattering coefficients %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpumpnormal = kx.*cos(al_grat) + kz.*sin(al_grat); % project pump propgation constant onto grating direction (in x,z plane)

alphas3 = sqrt(2*pi)./w0 ./ (2*cos(al_inc).^2).^2 .* (pi./Lambda).^2 .* (w./sin(2*theta)).^2 .* (dng./neff).^2;
alphas3 = alphas3 .* exp(-2*(w./sin(2*theta)).^2 .*(kpumpnormal.*cos(theta).^2-pi./Lambda).^2);
alphas3 = alphas3.*sin(phivar3d);

alphap3 = alphas3 .* cos(2*al_inc).^2;

if plotfig
    figure(26)
    plot(al_prop*180/pi,alphas3,al_prop*180/pi,alphap3,al_prop*180/pi,pxp.^2.*alphap3+pxs.^2.*alphas3,al_prop*180/pi,pyp.^2.*alphap3+pys.^2.*alphas3)
    xlabel('propagation angle')
    ylabel('scattering rate')
    legend('s pol','p pol','z pol, avg','y pol, avg')
    title('includes change in kx, incident angle, scatter power direction')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scattered electric field %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Es3 = -1 ./ (2*cos(al_inc).^2) .* w./sin(2*theta) .* dng./neff .* pi^1.5 ./Lambda;
Es3 = Es3 .* exp(-(w./(2*sin(2*theta))).^2 .*(2*kpumpnormal.*cos(theta).^2-2*pi./Lambda).^2);

Ep3 = -cos(2*al_inc) .* Es3 ;

% input vertical (y) polarisation

Eyx = Es3.*pys.*noutsx + Ep3.*pyp.*noutpx;
Eyy = Es3.*pys.*noutsy + Ep3.*pyp.*noutpy;
Eyz = Es3.*pys.*noutsz + Ep3.*pyp.*noutpz;

% input horizontal (x-z) polarisation

Exx = Es3.*pxs.*noutsx + Ep3.*pxp.*noutpx;
Exy = Es3.*pxs.*noutsy + Ep3.*pxp.*noutpy;
Exz = Es3.*pxs.*noutsz + Ep3.*pxp.*noutpz;

if plotfig
    figure(41)
    plot(al_prop*180/pi,Eyx,al_prop*180/pi,Eyy,al_prop*180/pi,Eyz)
    xlabel('propagation angle')
    ylabel('Output E-field in x,y,z')
    legend('x component','y component','z component')
    title('Output E-field for vertically polarised input')

    figure(42)
    plot(al_prop*180/pi,Exx,al_prop*180/pi,Exy,al_prop*180/pi,Exz)
    xlabel('propagation angle')
    ylabel('Output E-field in x,y,z')
    legend('x component','y component','z component')
    title('Output E-field for horizontally polarised input')
end

if pol==1   % pump polarisation horizontal
    al = pxp.^2.*alphap3+pxs.^2.*alphas3;
    Ex = Exx;
    Ey = Exz;
    Ez = Exy;
else        % poump polarization vertical
    al = pyp.^2.*alphap3+pys.^2.*alphas3;
    Ex = Eyx;
    Ey = Eyz;
    Ez = Eyy;
end

return

