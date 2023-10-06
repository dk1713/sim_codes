function [E_pump, E_target] = define_fields(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Default values
lam  = 780e-9; %#ok<NASGU> 
type = 'gaussian';
dire = [0, 1, 0];
view = 'horizontal';

switch nargin
    case 4
        lam=varargin{4};
    case 5
        lam=varargin{4}; type=varargin{5};
    case 6
        lam=varargin{4}; type=varargin{5}; dire=varargin{6};
    case 7
        lam=varargin{4}; type=varargin{5}; dire=varargin{6}; view=varargin{7};
    otherwise
        disp('Insufficient number of inputs!');
        return;
end

% propagation constant in air
n_eff   = 1.4635;
k1      = 2*pi*n_eff/lam;

% define the meshgrid
if strcmpi(view, 'horizontal') && length(varargin{2}) == 1
    [xx, zz] = meshgrid(varargin{1}, varargin{3});
    yy = varargin{2};
elseif strcmpi(view, 'vertical') && length(varargin{3}) == 1
    [xx, yy] = meshgrid(varargin{1}, varargin{2});
    zz = varargin{3};
else
    disp('Wrong viewpoint input or incorrect dimensions for y/z');
    return;
end

% define the pump
k1x = k1; w_y = 2e-6; w_z = 100*5.8e-6;
E_pump = exp(1i*k1x*xx).*exp(-zz.^2/w_z^2).*exp(-yy.^2/w_y^2);

% define the target field
n   = dire/norm(dire);
E0  = 0.1; % just amplitude.

if strcmpi(type, 'plane')
    E_target = E0*exp(1i*k1*(n(1)*xx + n(2)*yy + n(3)*zz));

elseif strcmpi(type, 'gaussian') %NB: working in the rotated frame.
    % beam waist
    w = 2.5e-6; %% ADJUST HERE %%
    % distance to target location
    zfoc = 50e-6; %% ADJUST HERE %%
    
    % z in rotated frame = projection onto prop. direction
    z = xx*n(1) + yy*n(2) + zz*n(3);
    % x=r in rotated frame = the orth. component to z
    r = sqrt(xx.^2 + yy.^2 + zz.^2 - z.^2);
    z = z - zfoc;
    
    zR = pi*w^2/lam;
    wz = w*sqrt(1+(z/zR).^2);
    eta = atan(z/zR);
    % inverse of curvature
    Rzi = z./(z.^2 + zR^2);
    
    E_target = E0 * w./wz .* exp(-(r./wz).^2) ...
        .* exp(1i*(k1*z+k1*r.^2.*Rzi/2-eta));
else
    disp('Not supported target field type.');
    return;
end

end