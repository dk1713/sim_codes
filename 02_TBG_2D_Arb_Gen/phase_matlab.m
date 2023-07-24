function out = phase_matlab(x_arg)
% A function handle saved in as .m file.
% global x phase
% 
% if (~exist('phase', 'var') || isempty(phase) )
%     load('phase.mat', 'x', 'phase');
% end
load('phase.mat', 'x', 'phase'); % global doesn't update!
out = interp1(x, phase, x_arg, 'spline', 0);
end