function r = cdfhat(y,h,x)
    k = (1/h)*(x-y);
    sumf =  0.75*((2/3)+k-(1/3)*k.^3).*(abs(k)<=1)+(k>1);
    r = mean(sumf);   
end

