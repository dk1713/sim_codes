function E = gaussian_on_surface(varargin)
%gaussian_on_surface: 
%   Detailed explanation goes here
% Assign input arguments
x = varargin{1}; lam = varargin{2}; n_prop = varargin{3}; 
pos_foc = varargin{4}; len_grat = varargin{5};
x_focus = pos_foc(1); y_focus = pos_foc(2);

if nargin == 5
    n = 1;
    is_focussed = 1;
elseif nargin == 6
    n = varargin{6};
    is_focussed = 1;
else
    n = varargin{6};
    is_focussed = varargin{7};
end

k_prop = 2*pi*n_prop /lam;

%% Focussed beam feature
% rotated frame
phi     = atan(x_focus/y_focus);
x_rot   = (x - x_focus)*cos(phi) - (0 - y_focus)*sin(phi);
y_rot   = (x - x_focus)*sin(phi) + (0 - y_focus)*cos(phi);

% Waist on the surface of the glass
w       = .5 * len_grat / exp(1);

% Rayleigh range(1)
temp_b  = pi*n_prop*w^2/lam;

if is_focussed == 1
    y_R     = .5*( temp_b - sqrt(temp_b^2 - 4*y_focus^2) );
else
    y_R     = .5*( temp_b + sqrt(temp_b^2 - 4*y_focus^2) );
end

% beam waist(1)
w_0     = real(sqrt( y_R * lam / (pi * n_prop) ));
w_rot   = w_0*sqrt(1 + (y_rot./y_R).^2);

% Radius of curvature(y)
R       = y_rot + y_R.^2./y_rot;

% Gouy phase
eta     = atan(y_rot./y_R);
q_inv   = 1./R - 1i*lam./(n_prop*pi*w_rot.^(2*n));

% Gaussian beam on the surface of the glass
E      = w_0./w_rot ...
    .*exp( -1i*(k_prop.*y_rot - eta) -1i*k_prop * x_rot.^(2*n) .* q_inv/2 );
end