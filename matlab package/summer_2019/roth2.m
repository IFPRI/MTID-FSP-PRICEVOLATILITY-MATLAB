function hrot = roth2(y,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    n=size(y,1);
    %aux = [ones(n,1) x 0.5*x.^2 (1/6)*x.^3 (1/24)*x.^4];
    aux = [ones(n,1) x 0.5*x.^2];
    b = (aux'*aux)\(aux'*y);
    s2=(1/(n-size(aux,2)))*(y-aux*b)'*(y-aux*b);
    C1K = (1/(2*sqrt(pi)))^(1/5);
    range = max(x)-min(x);
    theta22 = b(3)^2;
    %theta22 = mean(  (b(3)+b(4)*x+0.5*b(5)*x.^2).^2  );      
    hrot = n^(-0.2)*C1K*(range*s2/theta22)^(1/5);       
    
    

