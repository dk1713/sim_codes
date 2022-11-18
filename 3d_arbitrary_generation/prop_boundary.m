function [EE, k_cen_x, k_cen_y, phase_2] = prop_boundary(EE, xx, yy, k1, k2, k_x, k_y, kk_x, kk_y, n1, n2, dist, pol)
%UNTITLED Summary of this function goes here
%   Computing the propagated field by the path of 
sin_psi_x1  = k_x/k1;
sin_psi_y1  = k_y/k1;

sin_psi_x2  = n1*sin_psi_x1/n2;
sin_psi_y2  = n1*sin_psi_y1/n2;

k_cen_x     = k2*sin_psi_x2;
k_cen_y     = k2*sin_psi_y2;

% Computing for wave vectors and phases:
kk_t    = sqrt((kk_x + k_cen_x).^2 + (kk_y + k_cen_y).^2);

kk_z    = real( sqrt(k2^2 - kk_t.^2) );
phase_1 = asin(kk_t/(k1));
phase_2 = asin(kk_t/(k2));

% Uncomment the following to check if the central frequency is within the
% range acceptable.
% fprintf('Checking for limits of central wavenumber\n')
% fprintf('k_t(x) range: (%3.3e, %3.3e)',min(k_x+k_cen_x),max(k_x+k_cen_x))
% fprintf('\nk_cen_x  = %3.3e\n', k_cen_x)
% fprintf('k_t(y) range: (%3.3e, %3.3e)',min(k_y+k_cen_y),max(k_y+k_cen_y))
% fprintf('\nk_cen_y  = %3.3e\n', k_cen_y)
% fprintf('--------------------------------------------------------------\n')

EE    = EE .* exp(-1i*(k_cen_x*xx + k_cen_y*yy));
EE_k  = apply_Fresnel(fftshift(fft2(fftshift(EE))), ...
    phase_1, phase_2, n1, n2, pol)  ...
    .* exp(1i*kk_z.*dist); % propagation back onto the surface.
EE    = fftshift(ifft2(fftshift(EE_k)))...
    .* exp(1i*(k_cen_x*xx + k_cen_y*yy)); % add back central frequency.
end

