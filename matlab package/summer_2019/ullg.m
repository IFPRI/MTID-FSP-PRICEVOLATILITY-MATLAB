function m = ullg(y,x,h,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    n = size(y,1); 
    m = zeros(size(z,1),1); 
    b = gsubtract(x,z');
    a = (1/h)*b;
    ep = a.*(abs(a)<=1);
    k=0.75*(1-ep.^2);
    
    a1 = sum(k); 
    a2 = sum(k.*b); 
    a3 = sum(k.*(b.^2));
    b1 = sum(k.*y); b2 = sum((k.*b).*y);
    m = ((a3.*b1 - a2.*b2)./(a1.*a3-a2.^2+ n^(-2)))';
end

