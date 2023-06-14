function [alphas3, alphap3] = compute_alpha(varargin)
%   formula of power loss from tilted grating
%   input:  arrays (all of same size) of: 
%           local grating wavelength Lambda (m), 
%           grating direction algrat (rad), 
%           tilt angle (rad), 
%           transverse pump propagation constant kz (1/m)
%   one value for   pol=1 (input horizontal polarisation) or
%                   pol=2 (input vertical polarisation)
%   output: corresponding arrays:
%           al = energy loss rate (1/m)
%           Ex, Ey, Ez = scattered electrid field
% NB: in output x=pump propagation direction, y=in-plane transverse, 
% z=vertical but internally in this function y=vertical, z=transverse

    switch nargin
        case 9
        otherwise
            disp('Not enough input arguments');
            return;
    end
    
    Lambda  = varargin{1}; % grating period.
    al_grat = varargin{2}; % frame is rotated by this angle.
    theta   = varargin{3}; % tilt angle, angle of the slices to the y-axis
    w0      = varargin{4}; % waveguide width
    sigma   = varargin{5}; % input mode width
    beta    = varargin{6}; % propagation constant in the slab
    dng     = varargin{7}; % maximum dng taken for calculation
    neff    = varargin{8}; % effective index
    kz      = varargin{9}; % K_z

    % effective width of mode and grating
    w = 1/sqrt(1/w0^2+1/sigma^2);

    kx = sqrt(beta^2-kz.^2);       % longitudinal k
    al_prop = asin(kz/beta);       % propagation angle

    ngratx = -cos(theta).*cos(al_grat);   
    ngraty = sin(theta);  
    ngratz = -cos(theta).*sin(al_grat);   

    npropx = cos(al_prop);  npropy = 0*npropx;  npropz = sin(al_prop);

    al_inc = acos(-ngratx.*npropx -ngraty.*npropy -npropz.*ngratz);

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

    % input p pol is normal to propagation direction and s-pol vector, 
    % ninp ~ nins x nprop
    ninpx = ninsy.*npropz - ninsz.*npropy;
    ninpy = ninsz.*npropx - ninsx.*npropz;
    ninpz = ninsx.*npropy - ninsy.*npropx;

    % propagation direction of scattered beam 
    % nscat = nprop - 2*ngrat*(nprop.ngrat)
    qprojection = npropx.*ngratx + npropy.*ngraty + npropz.*ngratz;

    nscatx = npropx - 2*ngratx.*qprojection;
    nscaty = npropy - 2*ngraty.*qprojection;
    nscatz = npropz - 2*ngratz.*qprojection;

    % output s polarisation: same as input s polarisation
    noutsx = nincx; noutsy = nincy; noutsz = nincz;

    % p pol is normal to new propagation direction and s-pol vector,
    % noutp ~ nouts x nscat
    noutpx = noutsy.*nscatz - noutsz.*nscaty;
    noutpy = noutsz.*nscatx - noutsx.*nscatz;
    noutpz = noutsx.*nscaty - noutsy.*nscatx;

    % input vertical (y) polarisation
    % direction of polarisation
    nex = 0;
    ney = 1;
    nez = 0;

    % project of polarisation on s and p
    pys = nex.*ninsx + ney.*ninsy + nez.*ninsz;
    pyp = nex.*ninpx + ney.*ninpy + nez.*ninpz;

    % input horizontal (x-z) polarisation
    % ne = [-sin(al_prop), 0, cos(al_prop)];
    
    % direction of polarisation, orthogonal to nprop
    nex = -sin(al_prop);
    ney = 0;
    nez = cos(al_prop);

    % project of polarisation ons and p
    pxs = nex.*ninsx + ney.*ninsy + nez.*ninsz;
    pxp = nex.*ninpx + ney.*ninpy + nez.*ninpz;

    % output propagation direction and angle (from diffraction)
    Kgrat = 2*pi./Lambda;   % grating "propagation" constant

    kdiffx = kx - Kgrat.*cos(al_grat);      % diffracted beam
    kdiffz = kz - Kgrat.*sin(al_grat);
    kdiffy = sqrt(beta^2-kdiffx.^2-kdiffz.^2);

    phivar3d = real(asin(kdiffy/beta));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scattering coefficients %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % project pump propgation constant onto grating direction (in x,z plane)
    kpumpnormal = kx.*cos(al_grat) + kz.*sin(al_grat);
    
    alphas3 = sqrt(2*pi)./w0 ./ (2*cos(al_inc).^2).^2 ...
        .*(pi./Lambda).^2 .*(w./sin(2*theta)).^2 .*(dng./neff).^2 ...
        .* exp(-2*(w./sin(2*theta)).^2 ...
        .*(kpumpnormal.*cos(theta).^2 - pi./Lambda).^2);
    alphas3 = alphas3.*sin(phivar3d);
    
    alphap3 = alphas3 .* cos(2*al_inc).^2;
end