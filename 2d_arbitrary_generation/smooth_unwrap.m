function x = smooth_unwrap(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% check if the function is increasing or decreasing
if x(2) > x(1)
    is_inc = 1;
else
    is_inc = 0;
end

if is_inc
    for iter = 1:length(x) - 1
        while x(iter+1) <= x(iter)
            x(iter+1) = x(iter+1) + 2*pi;
        end 
    end
else
    for iter = 1:length(x) - 1
        while x(iter+1) >= x(iter)
            x(iter+1) = x(iter+1) - 2*pi;
        end 
    end
end
    
end

