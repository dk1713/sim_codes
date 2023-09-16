% generation of plane wave in arbitrary 3D direction with fixed polarisation
% by two tilted gratings for two pump beams (in x and y direction), both
% with horizontal polarisation

% based on grating_design3D_v1.m
% v2: allow for matrices as inputs (e.g. different angles)

global lam c k0 neff beta w0 sigma dng

% parameters
f
lam = 780e-9;           % wavelength
k0 = 2*pi/lam;

polout = 1;     % set output polarisation: 1=lin pol 1, 2=lin pol 2, 3=circ pol 1, 4=circ pol 2


c = 3e8;

% neff = 1;  %NOT IMPLEMENTED EVERYWHERE           % effective index of fundamental mode
neff = 1;
beta = neff*k0;         % propagation constant of mode
w0 = 2e-6;              % waist of fundamental mode (vertical)
sigma = 2e-6;           % waist of refractive index profile (vertical)
dng = 1.2e-3;           % grating index contrast (=1 for design calculation; realistic <=5e-3 ?)


% define target plane wave: direction, polarisation

%ntar = [1 1 3];   % direction of propagation (x=fwd, pump 1; y=transverse, pump 2; z=up)
%ntar = ntar/norm(ntar);

xx = -2:0.1:2;
yy = -1.9:0.1:1.9;
[ntarx,ntary] = meshgrid(xx,yy);
ntarz = 1 + 0*ntarx;
ntarnorm = sqrt(ntarx.^2+ntary.^2+ntarz.^2);
ntarx = ntarx./ntarnorm;
ntary = ntary./ntarnorm;
ntarz = ntarz./ntarnorm;


    % define two polarisation directions for output plane wave: 
    %   pol 1 orthogonal to propagation, in x,y plane
    %   pol 2 orhtogonal to propagation, orthogonal to pol 1

% if ntar(3)==1
%     pout1 = [0 -1 0]; 
%     pout2 = [1 0 0];
% else
%     pout1 = [ntar(2) -ntar(1) 0];    % ntar x [0 0 1]
%     pout1 = pout1/norm(pout1);
%     pout2 = [ntar(2)*pout1(3)-ntar(3)*pout1(2) ntar(3)*pout1(1)-ntar(1)*pout1(3)  ntar(1)*pout1(2)-ntar(2)*pout1(1)];  % ntar x ptar1
% end

pout1x = ntary;         % ntar x [0 0 1]
pout1y = -ntarx; 
pout1z = 0*ntarx;
pout1norm = sqrt(pout1x.^2+pout1y.^2+pout1z.^2);
pout1x = pout1x./pout1norm;
pout1y = pout1y./pout1norm;
pout1z = pout1z./pout1norm;

pout2x = ntarx.*pout1z-ntarz.*pout1y;    % ntar x ptar1
pout2y = ntarz.*pout1x-ntarx.*pout1z;  
pout2z = ntarx.*pout1y-ntary.*pout1x;

    % now correct for the upwards case (where the above is undefined)

pout1x(ntarz==1) = 0; 
pout1y(ntarz==1) = -1; 
pout1z(ntarz==1) = 0;
pout2x(ntarz==1) = 1; 
pout2y(ntarz==1) = 0;  
pout2z(ntarz==1) = 0;

figure(1); clf;
subplot(231)
pcolor(xx,yy,pout1x), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Pol 1 vector, x component')
subplot(232)
pcolor(xx,yy,pout1y), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Pol 1 vector, y component')
subplot(233)
pcolor(xx,yy,pout1z), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Pol 1 vector, z component')
subplot(234)
pcolor(xx,yy,pout2x), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Pol 2 vector, x component')
subplot(235)
pcolor(xx,yy,pout2y), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Pol 2 vector, y component')
subplot(236)
pcolor(xx,yy,pout2z), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Pol 2 vector, z component')


% set target polarisation: amplitudes (0,1) and relative phases (0,+/-pi/2)
% of the two output polarisation vectors

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

[lamgrat1,alphagrat1,alphatilt1] = grating_angles_3D_f2(ntarx,ntary,ntarz);
lamgrat1 = lamgrat1 * lam;
    
    % pump 2 (propagating in y): calculate in rotated frame (because grating_angles_3D_f2 assumes pump in x)

[lamgrat2,alphagrat2,alphatilt2] = grating_angles_3D_f2(ntary,-ntarx,ntarz);
lamgrat2 = lamgrat2 * lam;

figure(2); clf;
subplot(231)
pcolor(xx,yy,lamgrat1/lam), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Grating 1 period \Lambda/\lambda')
subplot(232)
pcolor(xx,yy,alphagrat1*180/pi), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Grating 1 direction')
subplot(233)
pcolor(xx,yy,alphatilt1*180/pi), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Grating 1 tilt')
subplot(234)
pcolor(xx,yy,lamgrat2/lam), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Grating 2 period \Lambda/\lambda')
subplot(235)
pcolor(xx,yy,alphagrat2*180/pi+90), shading flat, colorbar    % in real coordinate system
xlabel('x'), ylabel('y'), title('Grating 2 direction')
subplot(236)
pcolor(xx,yy,alphatilt2*180/pi), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Grating 2 tilt')


% calculate scattering rate for E-field and power (apart from dng factor, dng=1 here)

pol = 1;            % =1 horizontal input polarisation; =2 vertical input polarisation
ky = 0*lamgrat1;     % transverse prop. constant of pump

[al1,Ex1,Ey1,Ez1] = scatteringrate3_f(lamgrat1,alphagrat1,alphatilt1,ky,pol);
E1norm = sqrt(Ex1.^2+Ey1.^2+Ez1.^2);
%%
[al2,Ex2p,Ey2p,Ez2p] = scatteringrate3_f(lamgrat2,alphagrat2,alphatilt2,ky,pol); % in the wrong coordinate system
E2norm = sqrt(Ex2p.^2+Ey2p.^2+Ez2p.^2);

Ex2 = -Ey2p; Ey2 = Ex2p; Ez2 = Ez2p;    % scattered field from pump 2 in the true coordinate system


% calculate the required pump strengths (same dng) to achieve target polarisation

    % project scattered fields onto the two polarisation unit vectors
 
a1 = Ex1.*pout1x+Ey1.*pout1y+Ez1.*pout1z;    % scalar product E1.ptar1
a2 = Ex1.*pout2x+Ey1.*pout2y+Ez1.*pout2z;    % scalar product E1.ptar2

b1 = Ex2.*pout1x+Ey2.*pout1y+Ez2.*pout1z;    % scalar product E2.ptar1
b2 = Ex2.*pout2x+Ey2.*pout2y+Ez2.*pout2z;    % scalar product E2.ptar1

    % solve linear system for amplitudes of pumps 1 and 2

A1 = (b2.*ptar1-b1.*ptar2) ./ (a1.*b2-a2.*b1);
A2 = (-a2.*ptar1+a1.*ptar2) ./ (a1.*b2-a2.*b1);

figure(3); clf;   % remember: these pump strengths are for dng=1
subplot(221)
pcolor(xx,yy,abs(A1)), shading flat, colorbar
xlabel('x'), ylabel('y'), title(['Pol. ' num2str(polout) ': Pump 1 amplitude |A_1|'])
subplot(222)
pcolor(xx,yy,abs(A2)), shading flat, colorbar
xlabel('x'), ylabel('y'), title(['Pol. ' num2str(polout) ': Pump 2 amplitude |A_2|'])
subplot(223)
pcolor(xx,yy,angle(A1./A2)*180/pi), shading flat, colorbar
xlabel('x'), ylabel('y'), title('Relative phase of pumps')


% Crossterms: scattering of pump 1 on grating 2 and vice versa

    % scattering of pump 1 on grating 2

    % rotate grating 2 into coordinate system of pump 1
    
% if alphagrat2<0
%     alphagrat2p = alphagrat2 + pi/2;
%     upw1 = 1;    % =1 for pump 1 is scattered upwards by grating 2
% else
%     alphagrat2p = alphagrat2 - pi/2;
%     upw1 = 0;    % 0 for pump 1 is scattered downwards by grating 2
% end

upw1 = 1.0 * (alphagrat2<0);
alphagrat2p = alphagrat2 - pi/2 + pi*upw1;


    % calculate scattering rate and scattering direction
[al1cross,~,~,~] = scatteringrate3_f(lamgrat2,alphagrat2p,alphatilt2,ky,pol);

n1crossx = ( k0-2*pi./lamgrat2.*cos(alphagrat2p) ) / k0;
n1crossy = ( -2*pi./lamgrat2.*sin(alphagrat2p) ) / k0;
n1crossz = (-1).^(upw1+1).*sqrt(1-n1crossx.^2-n1crossy.^2);


    % scattering of pump 2 on grating 1

    % rotate grating 1 into coordinate system of pump 2
    
% if alphagrat1<0
%     alphagrat1p = alphagrat1 + pi/2;
%     upw2 = 0    % =1 for pump 2 is scattered upwards by grating 1
% else
%     alphagrat1p = alphagrat1 - pi/2;
%     upw2 = 1    % 0 for pump 2 is scattered downwards by grating 1
% end

upw2 = 1.0 * (alphagrat1>=0);
alphagrat1p = alphagrat1 + pi/2 - pi*upw2;


    % calculate scattering rate and scattering direction
[al2cross,~,~,~] = scatteringrate3_f(lamgrat1,alphagrat1p,alphatilt1,ky,pol);

n2crossy = ( k0-2*pi./lamgrat1.*cos(alphagrat1p) ) / k0;   % in pump 1 (true) coordinate system
n2crossx = ( +2*pi./lamgrat1.*sin(alphagrat1p) ) / k0;
n2crossz = (-1).^(upw2+1).*sqrt(1-n2crossx.^2-n2crossy.^2);


figure(4); clf;
subplot(231)
%pcolor(xx,yy,(real(n1crossx))), shading flat, colorbar
pcolor(xx,yy,n1crossx.*(abs(real(n1crossz))>0) ), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P1 on G2 scatter: direction, x component')
subplot(232)
%pcolor(xx,yy,(real(n1crossy))), shading flat, colorbar
pcolor(xx,yy,n1crossy.*(abs(real(n1crossz))>0) ), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P1 on G2 scatter: direction, y component')
subplot(233)
pcolor(xx,yy,(real(n1crossz))), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P1 on G2 scatter: direction, z component')
subplot(234)
%pcolor(xx,yy,(real(n2crossx))), shading flat, colorbar
pcolor(xx,yy,n2crossx.*(abs(real(n2crossz))>0) ), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P2 on G1 scatter: direction, x component')
subplot(235)
%pcolor(xx,yy,(real(n2crossy))), shading flat, colorbar
pcolor(xx,yy,n2crossy.*(abs(real(n2crossz))>0) ), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P2 on G1 scatter: direction, y component')
subplot(236)
pcolor(xx,yy,(real(n2crossz))), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P2 on G1 scatter: direction, z component')



figure(5); clf;
subplot(231)
pcolor(xx,yy,al1.*abs(A1).^2), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P1: design scatter \alpha|A|^2')
subplot(234)
pcolor(xx,yy,al1cross.*abs(A1).^2), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P1: parasitic scatter')
subplot(232)
pcolor(xx,yy,al2.*abs(A2).^2), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P2: design scatter \alpha|A|^2')
subplot(235)
pcolor(xx,yy,al2cross.*abs(A2).^2), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P2: parasitic scatter')
subplot(233)
pcolor(xx,yy,al1.*abs(A1).^2+al2.*abs(A2).^2), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P1+P2: design scatter \alpha|A|^2')
subplot(236)
pcolor(xx,yy,al1cross.*abs(A1).^2+al2cross.*abs(A2).^2), shading flat, colorbar
xlabel('x'), ylabel('y'), title('P1+P2: parasitic scatter')

