function [zhat,ez] = zhatq(p,n,o,l,q)
    zhat = o - (p(2)/p(1)).*( (n*(1-q)/l)^p(1) - 1 );
    ez = zhat*(1/(1+p(1)));  
end

