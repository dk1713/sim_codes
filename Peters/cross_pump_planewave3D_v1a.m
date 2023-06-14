% generation of plane wave in arbitrary 3D direction with fixed polarisation
% by two tilted gratings for two pump beams (in x and y direction), both
% with horizontal polarisation

% based on grating_design3D_v1.m
% v1a: allow for neff<>1 

global lam c k0 neff beta w0 sigma dng

% parameters

lam = 780e-9;           % wavelength
k0 = 2*pi/lam;

c = 3e8;

neff = 1;   % neff<>1 not really tested but "should" be ok now           % effective index of fundamental mode
beta = neff*k0;         % propagation constant of mode
w0 = 2e-6;              % waist of fundamental mode (vertical)
sigma = 2e-6;           % waist of refractive index profile (vertical)
dng = 1;                % grating index contrast (=1 for design calculation; realistic <=5e-3 ?)


% define target plane wave: direction, polarisation

ntar = [1 0 1];   % direction of propagation (x=fwd, pump 1; y=transverse, pump 2; z=up)
ntar = ntar/norm(ntar);

    % define two polarisation directions for output plane wave: 
    %   pol 1 orthogonal to propagation, in x,y plane
    %   pol 2 orhtogonal to propagation, orthogonal to pol 1

if ntar(3)==1
    pout1 = [0 -1 0]; 
    pout2 = [1 0 0];
else
    pout1 = [ntar(2) -ntar(1) 0];    % ntar x [0 0 1]
    pout1 = pout1/norm(pout1);
    pout2 = [ntar(2)*pout1(3)-ntar(3)*pout1(2) ntar(3)*pout1(1)-ntar(1)*pout1(3)  ntar(1)*pout1(2)-ntar(2)*pout1(1)];  % ntar x ptar1
end

polout = 1;     % set output polarisation: 1=lin pol 1, 2=lin pol 2, 3=circ pol 1, 4=circ pol 2

if polout==1    % 1=lin pol 1
    ptar1 = 1;
    ptar2 = 0;
end
if polout==2    % 2=lin pol 2
    ptar1 = 0;
    ptar2 = 1;
end
if polout==3    % 3=circ pol 1 
    ptar1 = 1;
    ptar2 = 1i;
end
if polout==4    % 4=circ pol 2
    ptar1 = 1;
    ptar2 = -1i;
end


% calculate grating properties (tilt, rotation, period) for the two pumps

    % pump 1 (propagating in x)

[lamgrat1,alphagrat1,alphatilt1] = grating_angles_3D_f2(ntar(1),ntar(2),ntar(3));
lamgrat1 = lamgrat1 * lam/neff;     % vXa: neff
    
    % pump 2 (propagating in y): calculate in rotated frame (because grating_angles_3D_f2 assumes pump in x)

[lamgrat2,alphagrat2,alphatilt2] = grating_angles_3D_f2(ntar(2),-ntar(1),ntar(3));
lamgrat2 = lamgrat2 * lam/neff;     % vXa: neff;
    
    
% calculate scattering rate for E-field and power (apart from dng factor, dng=1 here)

pol = 1;            % =1 horizontal input polarisation; =2 vertical input polarisation
ky = 0*lamgrat1;     % transverse prop. constant of pump

[al1,Ex1,Ey1,Ez1] = scatteringrate3_f(lamgrat1,alphagrat1,alphatilt1,ky,pol);
E1norm = sqrt(Ex1.^2+Ey1.^2+Ez1.^2);

[al2,Ex2p,Ey2p,Ez2p] = scatteringrate3_f(lamgrat2,alphagrat2,alphatilt2,ky,pol); % in the wrong coordinate system
E2norm = sqrt(Ex2p.^2+Ey2p.^2+Ez2p.^2);

Ex2 = -Ey2p; Ey2 = Ex2p; Ez2 = Ez2p;    % scattered field from pump 2 in the true coordinate system

% calculate the required pump strengths (same dng) to achieve target polarisation

    % project scattered fields onto the two polarisation unit vectors
 
a1 = Ex1*pout1(1)+Ey1*pout1(2)+Ez1*pout1(3);    % scalar product E1.ptar1
a2 = Ex1*pout2(1)+Ey1*pout2(2)+Ez1*pout2(3);    % scalar product E1.ptar2

b1 = Ex2*pout1(1)+Ey2*pout1(2)+Ez2*pout1(3);    % scalar product E2.ptar1
b2 = Ex2*pout2(1)+Ey2*pout2(2)+Ez2*pout2(3);    % scalar product E2.ptar1

    % solve linear system for amplitudes of pumps 1 and 2

A1 = (b2*ptar1-b1*ptar2) / (a1*b2-a2*b1)
A2 = (-a2*ptar1+a1*ptar2) / (a1*b2-a2*b1)

% Crossterms: scattering of pump 1 on grating 2 and vice versa

    % scattering of pump 1 on grating 2

    % rotate grating 2 into coordinate system of pump 1
    
if alphagrat2<0
    alphagrat2p = alphagrat2 + pi/2;
    upw1 = 1    % =1 for pump 1 is scattered upwards by grating 2
else
    alphagrat2p = alphagrat2 - pi/2;
    upw1 = 0    % 0 for pump 1 is scattered downwards by grating 2
end

    % calculate scattering rate and scattering direction
[al1cross,~,~,~] = scatteringrate3_f(lamgrat2,alphagrat2p,alphatilt2,ky,pol)

k1crossx = neff*k0-2*pi/lamgrat2*cos(alphagrat2p);     % vXa: neff
k1crossy = -2*pi/lamgrat2*sin(alphagrat2p);
k1crossz = (-1)^(upw1+1)*sqrt(neff^2*k0^2-k1crossx.^2-k1crossy.^2);     % vXa: neff
n1cross = [k1crossx k1crossy k1crossz]/(k0*neff)     % vXa: neff

    % scattering of pump 2 on grating 1

    % rotate grating 1 into coordinate system of pump 2
    
if alphagrat1<0
    alphagrat1p = alphagrat1 + pi/2;
    upw2 = 0    % =1 for pump 2 is scattered upwards by grating 1
else
    alphagrat1p = alphagrat1 - pi/2;
    upw2 = 1    % 0 for pump 2 is scattered downwards by grating 1
end

    % calculate scattering rate and scattering direction
[al2cross,~,~,~] = scatteringrate3_f(lamgrat1,alphagrat1p,alphatilt1,ky,pol)

% k2crossx = k0-2*pi/lamgrat1*cos(alphagrat1p);   % in pump 2 coordinate system
% k2crossy = -2*pi/lamgrat1*sin(alphagrat1p);
% k2crossz = (-1)^(upw2+1)*sqrt(k0^2-k2crossx.^2-k2crossy.^2);
% n2cross = [k2crossx k2crossy k2crossz]/k0

k2crossy = neff*k0-2*pi/lamgrat1*cos(alphagrat1p);   % in pump 1 (true) coordinate system     % vXa: neff
k2crossx = +2*pi/lamgrat1*sin(alphagrat1p);
k2crossz = (-1)^(upw2+1)*sqrt(neff^2*k0^2-k2crossx.^2-k2crossy.^2);     % vXa: neff
n2cross = [k2crossx k2crossy k2crossz]/(k0*neff)      % vXa: neff
