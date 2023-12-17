function EE_k = apply_Fresnel(EE_k, t1, t2, n1, n2, pol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if      pol == 's'
    EE_k = EE_k .* 2 * n1 .*cos(t1) ./ ( n1*cos(t1) + n2*cos(t2) );
elseif  pol == 'p'
    EE_k = EE_k .* 2 * n1 .*cos(t1) ./ ( n2*cos(t1) + n1*cos(t2) );
else
    fprintf('Error: incorrect output');
end

EE_k    = EE_k .* (imag(t1)==0);

end

