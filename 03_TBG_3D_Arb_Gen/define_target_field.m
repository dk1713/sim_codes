function E = define_target_field(varargin)
%   Defines the target field on the grating plane. Currently two cases
%   defined: case 7 = tilted Gaussian, and case 4/5 = plane wave
    switch nargin
        case 7 %case: tilted Gaussian
            type    = varargin{5};
            w0      = varargin{6}; % beam waist
            zfoc    = varargin{7}; % distance to focal point
        case 5 %case: plane wave
            type    = varargin{5};
        case 4 %case: plane wave (default)
            type = 'plane';
        otherwise
            disp('Not enough input arguments');
            return
    end
    xx      = varargin{1};
    yy      = varargin{2};
    k       = varargin{3};
    n_vec   = varargin{4};
    
    % target field amplitude - THIS SHOULD NOT HAVE ANY EFFECT
    E_0 = 0.1;      

    if strcmp(type, 'plane')   
        n_vec   = n_vec/norm(n_vec);
        k_x     = k*n_vec(1);
        k_y     = k*n_vec(2);    

        E = E_0*exp(1i*(k_x*xx+k_y*yy));
        
    elseif strcmp(type, 'gaussian')   
        n_vec = n_vec/norm(n_vec);

        % z in rotated frame = projection onto prop. direction
        zzr = xx*n_vec(1) + yy*n_vec(2);
        % x=r in rotated frame = the orth. component to z
        xxr = sqrt(xx.^2+yy.^2 - zzr.^2);
        zzr = zzr - zfoc;
        
        z_R = w0^2 * k/2;
        w_z = w0*sqrt(1+(zzr/z_R).^2);
        eta = atan(zzr/z_R);
        % inverse of curvature, 1/R_z
        Rz_inv = zzr./(zzr.^2 + z_R^2);

        E = E_0 * w0./w_z .* exp(-(xxr./w_z).^2) ...
            .* exp(1i*(k*zzr+k*xxr.^2.*Rz_inv/2-eta));
    else
        disp('Target field is not specified or incorrect input')
        return
    end
    
    % Check if the target field is correct. (uncomment below)
%     figure(101)
%     pcolor(x, y, abs(E).^2)    
%     xlabel('x'), ylabel('y')
%     title('target, |E|^2, z=0')
%     shading flat
%     axis equal
%     
%     figure(102)
%     pcolor(x, y, angle(E))    
%     xlabel('x'), ylabel('y')
%     title('target, phase(E), z=0')
%     shading flat
%     axis equal
end

