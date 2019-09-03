function [viol,p,color] = backtc(x,y,period,l,a)
% Inputs:
% x,y: are two column vectors of equal dimension
% period: is the size of subvectors of x and y whose elements will 
% be compared, elementwise.
% l: is the variable that indexes subvectors to be compared. l=1 starts
% with the first elements of x and y and compares all elements until
% position period.  l=2 starts with the second element and goes to position
% period+1, etc.
% a: a is a number between 0 and 1 to select the level of significance of
% the binomial test.
% Outputs: viol is the number of instances where the elements of x are
% greater than those of y.
% p is the p-value associated with the test
% 1-a is probability of "success", i.e., x>y 

ind=(x(l:period+l-1)>y(l:period+l-1));
viol=sum(ind);
p=2*(1-normcdf((abs(viol-period*(1-a)))/(sqrt(period*a*(1-a)))));

% If the p value is very small we are very confident this is unusual
if p<=.025
   color=2;
   
elseif .025<p && p<=.05
   color = 1;
   
else
   color=0;
end

