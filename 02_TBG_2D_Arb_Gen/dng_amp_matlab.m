function out = dng_amp_matlab(x_arg)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% global x dng_amp
% 
% if (~exist('dng_amp', 'var') || isempty(dng_amp) )
%     load('dng_amp.mat', 'x', 'dng_amp');
% end
load('dng_amp.mat', 'x', 'dng_amp');
out     = interp1(x, dng_amp, x_arg, 'spline', 0);
end