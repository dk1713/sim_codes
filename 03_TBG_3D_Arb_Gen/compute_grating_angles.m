function varargout = compute_grating_angles(varargin)
% input: 
%   n_outx          = direction of output field in x (1 or array),
%   n_outy          = direction of output field in y (1 or array),
%   n_outz          = direction of output field in z (1 or array),
%   n_in            = direction of input pump   (1x3)
%        where  x/y = pump propagation, y/x = transverse (in plane), 
%               z = vertical (out of plane)
% output: 
%   Lambda_grat     = grating period (units of lambda in material), 
%   alpha_ grat     = grating direction (rad), 
%   alpha_tilt      = tilt angle (rad)

    switch nargin
        case 4
            n_in = varargin{4};
        case 3
            n_in = [1, 0, 0];
        otherwise
            disp('Not enough input')
            return
    end

    % lengths normalised to wavelength (dimensionless)
    lam = 1;
    % k in material, assumed equal to beta
    k0 = 2*pi/lam;

    k_outn = sqrt(varargin{1}.^2 + varargin{2}.^2 + varargin{3}.^2);
    k_outx = k0*varargin{1}./k_outn;
    k_outy = k0*varargin{2}./k_outn;
    k_outz = k0*varargin{3}./k_outn;

    k_in    = k0*n_in;

    % grating K
    k_gratx = k_in(1) - k_outx;
    k_graty = k_in(2) - k_outy;

    % grating period, Lam_grat
    varargout{1} = 2*pi./sqrt(k_gratx.^2 + k_graty.^2);
    % grating direction (angle), alp_grat
    varargout{2} = atan2(k_graty, k_gratx);
%     varargout{2} = atan(k_graty./k_gratx);

    % normal to the grating planes
    n_gratx=k_outx-k_in(1); n_graty=k_outy-k_in(2); n_gratz=k_outz-k_in(3);
    n_gratz = n_gratz./sqrt(n_gratx.^2+n_graty.^2+n_gratz.^2);

    % grating tilt angle w.r.t. vertical, alp_tilt
    %   NB: (=0 vertical, =90 deg horizontal planes)
    varargout{3} = pi/2 - acos(n_gratz);
end
