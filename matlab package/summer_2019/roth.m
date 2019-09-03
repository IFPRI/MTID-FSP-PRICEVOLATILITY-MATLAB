function hrot = roth(y,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n=size(y,1);
    aux = [ones(n,1) x 0.5*x.^2 (1/6)*x.^3 (1/24)*x.^4];
    b = (aux'*aux)\(aux'*y);
    s2=(1/(n-size(aux,2)))*(y-aux*b)'*(y-aux*b);
    hrot = n^(-1/5)*( (s2*(max(x)-min(x))*( 1/(2*sqrt(pi)) ))/(mean((b(3)+b(4)*x+b(5)*0.5*x.^2).^2)) ).^(1/5);
end

