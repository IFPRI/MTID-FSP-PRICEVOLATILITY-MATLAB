clear all;
tic
yc=load('C:\projects\IFPRI\FSP\MATLAB\input\corn.txt');
N = size(yc,1);
disp(N)

lag2yc =yc(1:N-2);
lag1yc =yc(2:N-1);
yc=yc(3:N);

datac = [yc lag1yc lag2yc];
m = size(datac,1);

n=2000;
%repeatn=100;
repeatn=m-n;

store  = zeros(repeatn,11);
store1 = zeros(repeatn,7);
store2 = zeros(repeatn,7);
store3 = zeros(repeatn,7);
store4 = zeros(repeatn,7);

for repeati=1:repeatn
    y = datac(repeati:n+repeati-1,1);
    x = datac(repeati:n+repeati-1,2:3);
    evy = datac(n+repeati,1);
    evx = datac(n+repeati,2:3);
    repeati
    

    silvh = roth2(y,x(:,2));

    mhat  = mvllin(y,x,x,silvh);
    evmhat= mvllin(y,x,evx,silvh);
    
    silvhv = roth2((y-mhat).^2,x(:,2));
    
    shat2  = mvllin((y-mhat).^2,x,x,silvhv);
    evshat2= mvllin((y-mhat).^2,x,evx,silvhv);
    evshat2=0.00001*(evshat2<=0)+evshat2*(evshat2>0);
    
    yn=((y-mhat).*(shat2>0))./sqrt(shat2);

    
%     figure(1)
%     qqplot(yn)
%     grid on
%     box off
    
    l=round(n^(0.8-0.01));
    [xn,osn]=taildata(yn,l);
    
    afhatau = 1-(l/n);
    rothl = 0.79*n^(-0.2+0.01)*(iqr(yn));
    
    f = @(x)(cdfhat(yn,rothl,x)-afhatau);
    osns = fzero(f,osn);
    xns=yn-osns;
    xns=sortrows(xns((xns>0)));
   
    [phat,pci] = mle(xns,'pdf',@(xns,k,sigma,theta)gppdf(xns,k,sigma,0),'start',[0.15,.5],'Lowerbound',[0,0],'Upperbound',[0.25,1e25]);
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
    % store1 gives 95% quantile and CES
    store1(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(1,1)) (evmhat+sqrt(evshat2)*storeq(1,2)) ];
    % store2 gives 99% quantile and CES
    store2(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(2,1)) (evmhat+sqrt(evshat2)*storeq(2,2)) ];
    % store3 gives 99.5% quantile and CES
    store3(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(3,1)) (evmhat+sqrt(evshat2)*storeq(3,2)) ];
    % store4 gives 99.9% quantile and CES
    store4(repeati,:)=[ evy evx evmhat sqrt(evshat2) (evmhat+sqrt(evshat2)*storeq(4,1)) (evmhat+sqrt(evshat2)*storeq(4,2)) ];          
end
toc
period = 60;
corvec = zeros(repeatn-period+1,1);
indvec = zeros(repeatn-period+1,1);
pvec   = zeros(repeatn-period+1,1);

for l=1:repeatn-period+1
    [viol,p,color] = backtc(store1(:,1),store1(:,6),period,l,0.95);
    corvec(l) = color;
    indvec(l) = viol;
    pvec(l) = p;
end

t=(1:1:repeatn)';
xf1 = [t store1(:,1) store1(:,6) store1(:,7)];

% The graph 
figure(2)
plot(xf1(:,1),xf1(:,2),'black',xf1(:,1),xf1(:,3),'blue',xf1(:,1),xf1(:,4),'red')

% The next line saves figure(1) as a jpeg file called "plot", i.e.,
% "plog.jpeg" . Other formats are possible.

saveas(figure(2),'C:\projects\IFPRI\FSP\MATLAB\outputfolder\corn','png')

% The Output
% First column is a time stamp
% Second column is the time series of returns
% Third coliumns is VaR - 0.95
% Fourth Column is CES - 0.95
% Fifth Column is colora code

xf = [xf1(period:repeatn,1) xf1(period:repeatn,2) xf1(period:repeatn,3) xf1(period:repeatn,4) corvec];

% fprintf(fileID,'%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\n',xf');
% Next line saves xf in a ASCII file called "xf.txt"

%save('xf.txt','xf','-ascii')

fileID = fopen('C:\projects\IFPRI\FSP\MATLAB\outputfolder\output_corn.txt','w');
%fprintf(fileID,'%d\n',xf);
fprintf(fileID,'%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\n',xf');
fclose(fileID);

