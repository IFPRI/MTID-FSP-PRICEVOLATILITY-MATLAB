library cmlmt, pgraph;
#include cmlmt.sdf
graphset;


load retdata[] = D:\projects\IFPRI\FSP\gauss\input\corn.txt;
output file = D:\projects\IFPRI\FSP\gauss\outputfolder\output_corn.txt reset;
outwidth 256; /* to make sure that all outputs are printed in one line */
y=retdata[.,1];
"Number of returns" rows(y);


/*
lr1=lagn(r1,1);
lr2=lagn(r1,2);
dd_=r1~lr1~lr2;
dd_=dd_[3:rows(dd_),.];
*/

/*Number of observations used: window size*/
n=2000;
/*max rep=rows(y)-n-2.*/ 
rep=rows(y)-n-2; 
rows(y);

alph=0.05;

stor_z = zeros(rep,1);
stor_f = zeros(rep,1);
stor_xq = zeros(rep,1);
stor_xqn = zeros(rep,1);/*change*/

cycle2=1;
do while cycle2<=rep;
yt_2 = y[cycle2:n+cycle2-1];
yt_1 = y[cycle2+1:n+cycle2];
yt =y[cycle2+2:n+cycle2+1] ;
yf = y[n+cycle2+2];

q = { 0.05, 0.95 };
qdd_ = quantile(yt~yt_1~yt_2,q);
qdd_=qdd_';
reg = yt_1~yt_2;

/*ti = seqa(1,1,n);
xx = sortc(reg[.,1]~reg[.,2],1);
_pltype={6,6};
_plctrl = {-1,-1};
xy(xx[.,1],xx[.,2]);

_pltype = { 6 };
_plctrl = { 0 };
xy(ti,y);*/

ske = reg;

lbx= qdd_[2:3,1];
ubx= qdd_[2:3,2];

c = 0.4;
kn = minc(trunc(c*n^(2/5)*log(n))+1|trunc((0.5*n-1)/cols(reg)));
knots1 = seqa(lbx[1],(ubx[1]-lbx[1])/(kn+1),kn+2);
knots2 = seqa(lbx[2],(ubx[2]-lbx[2])/(kn+1),kn+2);

proc(1) = bspline(x, knots);
/* x is n by the number of conditioning variables in the E(Y|.) (d) and knots is (kn+2) by d. */
local n,d,z,i,j,kn;
n = rows(x); 
d = cols(x);
kn = rows(knots)-2;
z = zeros(n,d*kn);
i = 1; do while i<=d;
j = 1; do while j<=kn;
z[.,(j+ (i-1)*kn)] = ((x[.,i].>= (knots[j+1, i]*ones(n,1))) .and (x[.,i] .< (knots[j+2,i]*ones(n,1))));
j = j + 1;
endo;
i = i + 1;
endo;
retp(z);
endp;

xin = bspline(reg,(knots1~knots2));
xin1 = ones(n,1)~xin;
lambda = inv(xin1'*xin1)*xin1'*yt;

/*The next procedure provides estimates for the additive functions (equation (2.7) in
Wang and Yang).  The loop i is over all regressors (not including the intercept).  The output
of the procedure are d (nX1) columns with mhat_1,...,mhatd evaluated at the sample values of 
the regressors*/

proc(1) = dmean(x,l, kn);
local n,d, ma, mam,i;
n = rows(x); d = cols(x)/kn;
ma = zeros(n,d); mam = zeros(d,1);
i = 1;  do while i<=d;
mam[i] = meanc( x[.,((i-1)*kn+1):(i*kn)]*l[(1+(i-1)*kn+1):(1+i*kn)] );
ma[.,i] = xin[.,((i-1)*kn+1):(i*kn)]*l[(1+(i-1)*kn+1):(1+i*kn)] - mam[i];
i = i + 1;
endo;
retp(ma);
endp;

ma = dmean(xin,lambda,kn);

/*The following procedure calculates the estimated "pseudo-responses", i.e.,
\hat{y}_{id}.  The output is an n X d  matrix where the columns are 
\hat{y}_{i1},...,\hat{y}_{id}.  The procedure assumes d>=2*/

proc(1) = cest(y,m);
local n,d,yc,ind,i;
n = rows(y); d = cols(m);
yc = zeros(n,d); 
i = 1; do while i<=d;
if d == 1;
print "wrong inputs into dmean.";
endif;
if i == 1; 
ind = seqa(2,1,d-1);
elseif i == d;
ind = seqa(1,1,d-1);
else;
ind = seqa(1,1,i-1)|seqa(i+1,1,d-i);
endif;
yc[.,i] = y - meanc(y) - sumc((m[.,ind])');
i = i + 1;
endo;
retp(yc);
endp;

yia = cest(yt,ma);
hrotm = 4.25*stdc(reg)*n^(-1/5);

/*mah gives the estimated SBK-spline backfitted kernel estimators for 
each additive function.  In this code the second step is conducted 
using a Nadaraya-Watson estimator.  mat gives the oracle estimator*/

mah = zeros(n,cols(reg)); 

/*proc(1)= mnwgx(y,x,h,z);
local n,d,m_nw,mh,b,a,k,w,aux,j;
n=rows(y);d=cols(x);m_nw=zeros(rows(z),1); 
mh=diagrv(zeros(d,d),h); 
j=1; 
do while j<=rows(z);
b=(x'-ones(n,1)'.*.z[j,.]');
a=inv(mh)*b;
k=(abs(a).<=1).*(0.75*(ones(d,n)-a^2));
w=diagrv(zeros(n,n),prodc(k));
aux=ones(n,1);
m_nw[j]=inv(aux'*w*aux)*aux'*w*y;
j=j+1;
endo;
retp(m_nw);
endp;*/

proc nw(y,x,h,z);
local n,a,k,mhat;
    n=rows(y);
    a=(1/h)*(x-z');
    k=(abs(a).<=1).*(0.75*(ones(rows(x),rows(z))-a^2));
    mhat = sumc(k.*y)./sumc(k);
retp(mhat);
endp;

i = 1; 
do while i<=cols(reg);
mah[.,i] = nw(yia[.,i],reg[.,i],hrotm[i],reg[.,i]);
i = i + 1;
endo;

/*Graphs for the estimated functions*/

/*gr1 = reg[.,1]~mah[.,1]~(yt-meanc(yt)-mah[.,2]);
gr1 = sortc(gr1,1);
_pltype ={ 6,6 };
_plctrl ={ 0,-1};
xy(gr1[.,1],gr1[.,2]);

gr2 = reg[.,2]~mah[.,2]~(yt-meanc(yt)-mah[.,1]);
gr2 = sortc(gr2,1);
_pltype ={ 6,6 };
_plctrl ={ 0,-1};
xy(gr2[.,1],gr2[.,2]);*/

kn = minc((trunc(c*n^(2/5)*log(n))+1)|trunc((n/2-1)/cols(ske)));

xins = bspline(ske,(knots1~knots2));
xins1 = ones(n,1)~xins;
yst2 = (yt - meanc(yt)- sumc(mah'))^2;

lambdas = inv(xins1'*xins1)*xins1'*yst2;
sa = dmean(xins,lambdas,kn);
yias = cest(yst2,sa);

hrots = 5.25*stdc(ske)*n^(-1/5);
sah = zeros(n,cols(ske)); 

i = 1; 
do while i<=cols(ske);
sah[.,i] = nw(yias[.,i],ske[.,i],hrots[i],ske[.,i]);
i = i + 1;
endo;

/*Graphs for the estimated functions*/

/*gs1 = ske[.,1]~sah[.,1]~(yst2-meanc(yst2)-sah[.,2]);
gs1 = sortc(gs1,1);
_pltype={6};
_plctrl = {0};
xy(gs1[.,1],gs1[.,2]);

gs2 = ske[.,2]~sah[.,2]~(yst2-meanc(yst2)-sah[.,1]);
gs2 = sortc(gs2,1);
_pltype={ 6 };
_plctrl = { 0 };
xy(gs2[.,1],gs2[.,2]);*/

/* Estimated Variance */
evar = meanc(yst2)+sumc(sah');
evar = 0.0001*(evar .<=0)+evar.*(evar .>0);

xn = (yt - meanc(yt)- sumc(mah'))./sqrt(evar);

/*taildata(x,l) is a data procedure.  It gets the data on the residuals and produces a
new data set based of a fixed number of exceedances, given by "l", as well as
the lth order statistic, which is the random threshold for the tail*/

proc(2) = taildata(x,l);
local n,sdata,os,y;
n=rows(x);
sdata = sortc(x,1);
os = sdata[n-l,1];
y = sdata[n-l+1:n,1]-os;
retp(y,os);
endp;

/* l is the number of exceedances. It can be set at different levels.A condition is that 
l>n(1-q)=n*alph, where q is the quantile (1-alph)*/


	proc(2) = taildatan(x,l);/*change*/
		local n,sdata,os,y;/*change*/
		n=rows(x);/*change*/
		sdata = sortc(x,1);/*change*/
		os = sdata[l+1,1];/*change*/
		y = os-sdata[1:l,1];/*change*/
		retp(y,os);/*change*/
	endp;	/*change*/

l=200;

{ xn_n,osn_n } = taildatan(xn,l);/*change*/
{ xn,osn } = taildata(xn,l);     /*change*/

/*lmomgpd(x) is a procedure to calculate L-Moments estimators*/

proc lmomgpd(x);
local n,ct,sx,b0,b1,b2,l1,l2,k,sigma,par;
n=rows(x);
ct = seqa(1,1,n);
sx = sortc(x,1);
b0 = (1/n)*sumc(sx);
b1 = (1/(n*(n-1)))*sumc(ct[1:n-1].*sx[2:n]);
b2 = (1/(n*(n-1)*(n-2)))*sumc(ct[2:n-1].*ct[1:n-2].*sx[3:n]);
l1=b0; l2=2*b1-b0;
k=l1/l2 - 2;
sigma = (1+k)*l1;
par = sigma|k;
retp(par);
endp;

parn = lmomgpd(xn);
parn_n = lmomgpd(xn_n);/*change*/
/* Maximum Likelihood Estimation */

/* Data structure */
struct DS d0;
d0 = reshape(dsCreate,1,1);
d0[1].dataMatrix = xn;

/* Parameter structure */
struct PV p1;
p1 = pvCreate;
p1 = pvPacki(p1,parn,"pvar",1);

/* cmlmtControl Instance. */
struct cmlmtControl c1;
c1 = cmlmtControlCreate;
c1.Title = "Parameter estimates with MLE";
c1.Bounds = {0.001 1e256, 
             -1e256 1e256};

/* output for CML procedure. */
struct cmlmtResults out1;
/*out1 = cmlmt(&llgpd,p1,d0,c1); */
bini = pvUnpack(out1.par,1);
phn = bini;

proc llgpd(struct PV p1,struct DS d0,ind);
local sigma,k,z,e,logf,x,pvar,w,g;
struct modelResults logl;
x = d0[1].dataMatrix;
pvar = pvUnpack(p1,1);
sigma=pvar[1]; k=pvar[2];
z = x./ sigma;
if k .> 0;
e = (0 .<= x .and x .<= sigma ./ k);
    logf = - ln(sigma)- (1-k).*(-(1./k).*ln(1- k*z));
    logf= logf.*e + 0 .*(1-e);
elseif k .< 0;
e = (0 .< x) ;
    logf = - ln(sigma)- (1-k).*(-(1./k).*ln(1- k*z));
    logf= logf.*e + 0 .*(1-e);
    else;
logf = - ln(sigma)- z;
endif;
w = 1- k.*z;
w = missex(w,w.<= 0);
if k==0;
g = (1 ./ sigma).*(z-1);
else;
g = ones(rows(x),2);
g[.,1] = (x.*(1-k)./(sigma.*w) - 1)./sigma;
g[.,2] = ( z./(k.*w) + ln(w)./(k^2) ).*(k-1)-(1 ./k).*ln(w);
g = missrv(g,1);
endif;
if ind[1];
logl.Function = logf;
endif;
if ind[2];
logl.Gradient = g;
endif;
retp(logl);
endp;


/*VaR and Expected Shortfall Calculation*/

proc(2) = zhatq(par,n,os,q);
local zhat,ez;
zhat = os - (par[1]/par[2]).*( (n*(1-q)/l)^par[2] - 1 );
ez = zhat*(1/(1+par[2])+(par[1]+par[2]*os)./((1+par[2])*zhat));
retp(zhat,ez);
endp;

qalph = 1-alph;

{ lzhatn, lezn } = zhatq(parn,n,osn,qalph);
{ lzhatn_n, lezn_n } = zhatq(parn_n,n,osn_n,qalph);/*change*/
/*{ mzhatn, mezn } = zhatq(phn,n,osn,alphaq);*/

mf = zeros(1,cols(reg));
regf = yt[n]~yt_1[n];

i = 1; 
do while i<=cols(reg);
mf[.,i] = nw(yia[.,i],reg[.,i],hrotm[i],regf[.,i]);
i = i + 1;
endo;

sf = zeros(1,cols(reg));

i = 1; 
do while i<=cols(ske);
sf[.,i] = nw(yias[.,i],ske[.,i],hrots[i],regf[.,i]);
i = i + 1;
endo;


/* Estimated Variance */
evarf = meanc(yst2)+sumc(sf[1,.]');
evarf = 0.0001*(evarf .<=0)+evarf.*(evarf.>0);

xqnl =   meanc(yt) + sumc(mf[1,.]') + sqrt(evarf)*lzhatn;
xqnl_n =   meanc(yt) + sumc(mf[1,.]') + sqrt(evarf)*lzhatn_n;/*change*/

 stor_z[cycle2]= lzhatn;
stor_xq[cycle2]= xqnl;
stor_xqn[cycle2]= xqnl_n;/*change*/
stor_f[cycle2]=yf; 

cycle2 = cycle2+1;
endo;

yf=stor_f[1:rep];
xqnl=stor_xq[1:rep];
lzatn=stor_z[1:rep];

/*Backtest.g: x is the return series, and y is the quantile, a=alph*/

{indxqnl,pvalxqnl,inxqnl}=backt(yf,xqnl,alph);

proc(3)= backt(x,y,a);
    local n,m,ind,in,p;
    n=rows(x);
    m=rows(y);
    in=zeros(n,1);
    ind=(x[1:n].>y[1:m]);
    in[1:n,1]=ind;
    ind=sumc(ind);
    @for two-sided binomial test, a level=0.95.Normal approximation
    is used for the binomial distribution.Compare the p-value with 5%.@
    p=2*(1-cdfn( (abs(ind-n*a))/sqrt(n*a*(1-a)) ));
    retp(ind,p,in);
endp;

print "Backtest for n="rows(yf) "Alpha=" alph; 
print "The expected violation number is n times alpha =" rows(yf)*alph;

tres1=lzhatn';
tres2=ones(rows(tres1),1)*(rows(yf)*alph);
tres3= (indxqnl)';
tres4=(pvalxqnl)';
print "The residual quantile is" tres1;
print "The actual number of forecasted violations is " tres3;
print "The p-value is" tres4;


/*print tres1~tres2~tres3~tres4;*/
struct plotcontrol myplot;
myplot = plotGetDefaults("xy");


t=seqa(1,1,rows(yf));
plotSetXLabel(&myplot,"time");
plotSetYLabel(&myplot,"log returns");

// xlabel("time");
// ylabel("log returns");

_pltype = { 6 };
gr=sortc(t~yf~stor_xq[.,1]~stor_xqn[.,1],1);/*change*/
//xy(gr[.,1],gr[.,2]~gr[.,3]);
plotxy(myplot, gr[.,1],gr[.,2]~gr[.,3]~gr[.,4]);
plotSave("D:\\projects\\IFPRI\\FSP\\gauss\\outputfolder\\corn.jpg", 30|18);
//graphprt("-c=5 -cf=D:\\projects\\IFPRI\\FSP\\gauss\\outputfolder\\rice.bmp");

_ptek = "graph_corn.tkf";

/*print gr;*/

proc(4)= backtc(x,y,period,l,a);
    local n,m,ind,in,p,cor;
    n=rows(x);
    m=rows(y);
    in=zeros(n,1);
    ind=(x[l:period+l-1].>y[l:period+l-1]);
    ind=sumc(ind);
    @for two-sided binomial test, a level=0.95.Normal approximation
    is used for the binomial distribution.Compare the p-value with 5%.@
    p=2*(1-cdfn( (abs(ind-period*(1-a)))/sqrt(period *a*(1-a)) ));
    if p.<=0.025;/*The null that violations are consistent with expected violations is highly questionable (red)*/
    cor =2;
    elseif 0.05 .>= p .>= 0.025; /*the null that violations are consistent with expectations is questionable at a low level (orange)*/
    cor=1;
    else;
    cor = 0; /*accept the null that violations are consistent with expectations (green)*/
    endif;
    retp(ind,p,in,cor);
endp;
n=rows(gr);
period=60; corvec=zeros(n-period+1,1); indvec=zeros(n-period+1,1); pvec=zeros(n-period+1,1);

l=1;
do while l<=n-period+1;
{indr,pvalr,inr,corx}=backtc(gr[.,2],gr[.,3],period,l,0.95);
corvec[l]=corx;
indvec[l]=indr;
pvec[l]=pvalr;
l=l+1;
endo;

alpha=0.05;


vecviol = ones(n-period+1,1)*period*alpha;
print "First column is a time index,
second column are returns,
third column is the 95 percent conditional quantile
fourth column is the 5 percent conditional quantile,/*change*/
fifth column is the color code";
print gr[period:n,1]~gr[period:n,2]~gr[period:n,3]~gr[period:n,4]~corvec;/*change*/
