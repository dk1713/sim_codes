function dy = c_diff(x, y)
% Calculating the central difference first order derivatives with forward
% and backward difference for the boundary. NOTE: The end derivatives are
% not necessary for the overall calculation.

dy = zeros(size(y));
dx = x(2) - x(1);

dy(1) = (y(2) - y(1))/dx;

for iter = 2:length(x)-1
    dy(iter) = .5* (y(iter+1) - y(iter-1)) / dx;
end

dy(end) = (y(end) - y(end-1))/dx;

end

