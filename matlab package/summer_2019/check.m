clear all;
X =randn(2000,2);
e = randn(2000,1);
Y= 2*X(:,1).^2+X(:,2)+e;


silvh = roth2(Y,X(:,2));
mhat= mvllin(Y,X,X,silvh);

mhat