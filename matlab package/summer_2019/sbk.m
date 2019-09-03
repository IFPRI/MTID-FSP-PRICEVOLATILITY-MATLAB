clear all;
yc=load('corn.txt');
N = size(yc,1);

lag1yc = yc(1:N-1);
yc = yc(2:N);

datac = [yc lag1yc];
m = size(datac,1);

n=2000;
%repeatn=1;
repeatn=m-n;

store=zeros(repeatn,11);
store1 = zeros(repeatn,6);
store2 = zeros(repeatn,6);
store3 = zeros(repeatn,6);
store4 = zeros(repeatn,6);

for repeati=1:repeatn
    y = datac(repeati:n+repeati-1,1);
    x = datac(repeati:n+repeati-1,2);
    evy = datac(n+repeati,1);
    evx = datac(n+repeati,2);
    repeati
    silvh = roth(y,x);
    mhat=ullin(y,x,x,silvh);
    evmhat=ullin(y,x,evx,silvh);
    
    silvhv = roth((y-mhat).^2,x);
    shat2=ullin((y-mhat).^2,x,x,silvhv);
    evshat2=ullin((y-mhat).^2,x,evx,silvhv);
    evshat2=0.00001*(evshat2<=0)+evshat2*(evshat2>0);
    
    yn=((y-mhat).*(shat2>0))./sqrt(shat2);
    l=round(n^(0.8-0.01));
    [xn,osn]=taildata(yn,l);
    
    afhatau = 1-(l/n);
    rothl = 0.79*n^(-0.2+0.01)*(iqr(yn));
    
    f = @(x)(cdfhat(yn,rothl,x)-afhatau);
    osns = fzero(f,osn);
    xns=yn-osns;
    xns=sortrows(xns((xns>0)));
    
    [phat,pci] = mle(xns,'pdf',@(xns,k,sigma,theta)gppdf(xns,k,sigma,0),'start',[0.15,1]);
    k=-phat(1);
    sigma=phat(2);
    par = [k;sigma];
    
    alphavalue = [0.05;0.01;0.005;0.001];
    storeq=zeros(4,2);
    
    for alphai=1:size(alphavalue,1) 
        q=1-alphavalue(alphai);
        [zhat,ez]=zhatq(par,n,osns,l,q);
        storeq(alphai,:)=[zhat,ez];
    end
    
    store(repeati,:)=[ par' osns storeq(:,1)' storeq(:,2)' ];
    store1(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(1,1)) (evmhat+sqrt(evshat2)*storeq(1,2)) ];
    store2(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(2,1)) (evmhat+sqrt(evshat2)*storeq(2,2)) ];
    store3(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(3,1)) (evmhat+sqrt(evshat2)*storeq(3,2)) ];
    store4(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(4,1)) (evmhat+sqrt(evshat2)*storeq(4,2)) ];
          
end

period = 60;
corvec = zeros(repeatn-period+1,1);
indvec = zeros(repeatn-period+1,1);
pvec = zeros(repeatn-period+1,1);

for l=1:repeatn-period+1
    [viol,p,color] = backtc(store1(:,1),store1(:,5),period,l,0.95);
    corvec(l) = color;
    indvec(l) = viol;
    pvec(l) = p;
end

save_var = [indvec pvec corvec];
save('save_var')

% t=(1:1:repeatn)';
% xf1 = [t store1(:,1) store1(:,5) store1(:,6)];
% 
% figure(1)
% plot(xf1(:,1),xf1(:,2),'black',xf1(:,1),xf1(:,3),'blue',xf1(:,1),xf1(:,4),'red')

