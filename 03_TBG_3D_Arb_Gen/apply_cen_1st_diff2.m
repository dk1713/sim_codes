function f_out = apply_cen_1st_diff2(varargin)
%   Applying the first order differentiation using the central finite
%   difference approach in 2D matrix. col_or_row determines whether
%   difference is performed in column-wise or row-wise.
    switch nargin
        case 3 %case: tilted Gaussian
            f_in        = varargin{1};
            h           = varargin{2};
            col_or_row  = varargin{3};
        otherwise
            disp('Not enough input arguments');
            return
    end
    
    if col_or_row == 1
        f_out = [...
            -3*f_in(:,1) + 4*f_in(:,2) - f_in(:,3) ...
            f_in(:,3:end) - f_in(:,1:end-2) ...
            f_in(:,end-2) - 4*f_in(:,end-1) + 3*f_in(:,end)] / (2*h);
    elseif col_or_row == 2
        f_out = [...
            -3*f_in(1,:) + 4*f_in(2,:) - f_in(3,:);...
            f_in(3:end,:) - f_in(1:end-2,:);...
            f_in(end-2,:) - 4*f_in(end-1,:) + 3*f_in(end,:)] / (2*h);
    else
        disp('row_or_col is not inputted correctly!')
        return
    end
end

