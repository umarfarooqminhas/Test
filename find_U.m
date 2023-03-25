function [k] = find_U (xv,x)

x_max = max(xv);

if x_max >= x

    for i = 1 : length(xv)
        Dx(i) = abs (x- xv(i));
    end

    Dx_min = min(Dx);

    k = find(Dx_min == Dx);

else
    error('out of range')
end

 