% This function performs multivariate local linear estimation
% based an Epanechnikov kernel.
%INPUTS:
% y: an (nX1) column vector containing the observations on the regressand
% x: an (nXnreg) matrix containing the observations on the regressors
% evl: a (pXnreg) matrix containing the points where it is desired to
% evaluate the estimator.  p can be any integer {1,2,...}
% h: is a (nregX1) vector of bandwidths.

function mvllin = mvllin(y,x,evl,h)
n=size(y,1);
p=size(evl,1);
nreg=size(x,2);
K=ones(n,p);
s=zeros(p,n);
e = [1 zeros(1,nreg)];
for j=1:nreg
    a=(1/h)*gsubtract(x(:,j),evl(:,j)');
    %k=0.75*(1-a.^2).*(abs(a)<=1);
    k=(1/sqrt(2*pi))*exp(-0.5*a.^2);
    K=K.*k;
end    

for i=1:p
    b=gsubtract(x,evl(i,:));
    R=[ones(n,1) b];
    %P=diag(K(:,i));
    Q=R.*K(:,i);
    %s(i,:)=e*((R'*P*R)\(R'*P));
    s(i,:)=e*((Q'*R)\(Q'));
end
mvllin = s*y;