% This function performs multivariate local linear estimation
% based an Epanechnikov kernel.
%INPUTS:
% y: an (nX1) column vector containing the observations on the regressand
% x: an (nXnreg) matrix containing the observations on the regressors
% evl: a (pXnreg) matrix containing the points where it is desired to
% evaluate the estimator.  p can be any integer {1,2,...}
% h: is a (nregX1) vector of bandwidths.

function m = ullin(y,x,evl,h)
n=size(y,1);
b = gsubtract(x,evl');
a = (1/h)*gsubtract(x,evl');
k=(1/sqrt(2*pi))*exp(-0.5*a.^2);
%k = (abs(a) <= ones(n,size(a,2))).*(0.75*(ones(n,size(a,2))-a^2));
a1 = sum(k); 
a2 = sum(k.*b); 
a3 = sum(k.*b.^2);
b1 = sum(k.*y); 
b2 = sum((k.*b).*y);
m = ((a3.*b1 - a2.*b2)./(a1.*a3-a2.^2+ n^(-2)))';



 








