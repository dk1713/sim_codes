function EE_k = apply_Fresnel(EE_k, varargin)
%apply_Fresnel: Takes in electric field of any dimensions and apply
%Fresnel's coefficients to the transmission part of the propagation. The
%propagation is transversing from medium 1 to 2.
%   Input:
%   EE_k    (nxm float) electric field in Fourier domain in medium 1.
%   t1      (nxm float) phase of medium 1.
%   t2      (nxm float) phase of medium 2.
%   n1      (float) index of medium 1.
%   n2      (float) index of medium 2.
%   pol     (string) indication for s- or p-polarisation.
%   Output:
%   EE_k    (nxm float) electric field in Fourier domain in medium 2.

    switch nargin
        case 5
            pol = 's'; % default pol
        case 6 
            pol = varargin{5};
        otherwise
            disp('Error: not enough arguments');
            return
    end
    
    t1 = varargin{1}; t2 = varargin{2}; 
    n1 = varargin{3}; n2 = varargin{4};

    if      pol == 's'
        EE_k = EE_k .* 2 * n1 .*cos(t1) ./ ( n1*cos(t1) + n2*cos(t2) );
    elseif  pol == 'p'
        EE_k = EE_k .* 2 * n1 .*cos(t1) ./ ( n2*cos(t1) + n1*cos(t2) );
    else
        fprintf('Error: invalid polarisation\n');
        return
    end

    EE_k    = EE_k .* (imag(t1)==0);
    return
end

