function [half_max, x_at_hm, fwhm] = find_fwhm(x_arg, y_arg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
half_max    = .5*max(y_arg);

iter        = 1;
x_at_hm     = zeros(1,2);
counter     = 1;

while  iter < length(x_arg) && counter < 3
    if 100*abs(half_max - y_arg(iter))/half_max < 1
        x_at_hm(counter) = x_arg(iter);
        counter = counter + 1;
        iter = iter + floor(.1*length(x_arg));
    end
    iter = iter + 1;
end

fwhm        = x_at_hm(2) - x_at_hm(1);

end

