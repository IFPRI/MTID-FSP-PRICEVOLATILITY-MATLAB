function [y,os] = taildata(x,l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    n=size(x,1);
    sdata = sortrows(x);
    os = sdata(n-l);
    y = sdata(n-l+1:n)-os;
end

