function iqr = iqr(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    n = size(x,1);
    s = sortrows(x);
    iqr = s(round(0.75*n))-s(round(0.25*n));
end

