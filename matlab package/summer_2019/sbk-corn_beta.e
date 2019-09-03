/*Call relevant libraries*/

library pgraph,cmlmt,nlsysmt,co;
#include cmlmt.sdf
#include nlsysmt.sdf;
#include co.ext;
coset;
graphset; 
struct nlControl nlc;
struct nlOut nlo;
nlc = nlControlCreate;

/* Next two lines loads return data */
filestr = "corn10272014.txt";
load retdata[] = D:\projects\IFPRI\FSP\gauss\input\corn.txt;

/* Next two names the output file */
filestr = "outcorn10272014_empirical.txt";
output file = D:\projects\IFPRI\FSP\gauss\outputfolder\output_corn_beta.txt reset;

yc = retdata[.,1];
datac = yc~lag(yc);
datac = datac[2:rows(datac),.];

/* m is the total number of returns.  All historic series. */
m = rows(datac);

/* n is the sample size to be used in the estimation */
n=2000;

repeatn = m-n; 
/*repeatn=100;*/

store = zeros(repeatn, 11);
storey1 = zeros(repeatn,6);
storey2 = zeros(repeatn,6);
storey3 = zeros(repeatn,6);
storey4 = zeros(repeatn,6);

repeati = 1; 
do while repeati<=repeatn;
    
    y = datac[repeati:(n+repeati-1),1]; 
    x = datac[repeati:(n+repeati-1),2]; 
    evy = datac[n+repeati,1]; 
    evx = datac[n+repeati,2];


/*Bandwith ROT Selection procedure*/ 
   proc roth(y,x);
    local hrot,n,aux,b,s2;
            n=rows(y);
            aux = ones(n,1)~x~0.5*x^2~inv(6)*x^3~inv(24)*x^4;
            b = inv(aux'*aux)*aux'*y;
            s2=inv(n-cols(aux))*(y-aux*b)'*(y-aux*b);
            /* hrot = n^(-0.2)*(s2*(maxc(x)-minc(x))*( (2*sqrt(pi))^(-1) )*inv( n^(-1)*sumc((b[3]+b[4]*x+b[5]*0.5*x^2)^2)))^(1/5);*/
            hrot = n^(-0.2)*(s2*(maxc(x)-minc(x))*( 15 )*inv( n^(-1)*sumc((b[3]+b[4]*x+b[5]*0.5*x^2)^2)))^(1/5);
    retp(hrot);
    endp;

/* Next two lines calculate bandwidth using roth and estimates regression at all data points. */
    silvh = roth(y,x); 
    mhat = ullg(y,x,silvh,x);

/* Next two lines calculate bandwidth using roth and estimates the skedastic function at all data points. */
    silvhv = roth((y-mhat)^2,x);
    shat2 = ullg((y-mhat)^2,x,silvhv,x);

/* Next two lines estimate the regression and the skedastic function at evx (one time epriod ahead). */
    evmhat = ullg(y,x,silvh,evx);
    evshat2 = ullg((y-mhat)^2,x,silvhv,evx);

/* The next procedure accounts for the possibility of the variance being negative at the point of evaluation. */
    evshat2 = 0.0001*(evshat2.<=0)+evshat2.*(evshat2.>0);

/* The next procedure accounts for the possibility of the variance being negative at the sample points. */
    yn = (y - mhat);
    if minc(shat2).<=1e-4;
        indzsv = {-1e256,1e-4};
        indzs = indexcat(shat2,indzsv);
        yn = putvals(yn,indzs, zeros(rows(indzs),1));
        shat2p = putvals(shat2,indzs, ones(rows(indzs),1));
        yn = yn./sqrt(shat2p);
    elseif minc(shat2).>1e-4;
        yn = yn./sqrt(shat2);
    endif;
/* yn are the standardized residulas of equation (4) in Martins-Filho et al (2018)*/
   
/* The next two procedure are not used.*/
proc cdfcauchyi(p,m,s);
    local u;
    u = m + s*tan(pi*(p-0.5));
    retp(u);
endp;
/*lmomgpd(x) is a procedure to calculate L-Moments estimators. In all "x" is data to be used, it is the result of the taildata procedure.*/
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

/* taildata(x,l) is a data procedure. It gets the data on the residuals and produces a new data set based of a fixed number of exceedances, 
given by "l", as well as the lth order statistic, which is the random threshold for the tail.*/

proc(2) = taildata(x,l);
    local n,sdata,os,y;
    n=rows(x);
    sdata = sortc(x,1);
    os = sdata[n-l,1];
    y = sdata[n-l+1:n,1]-os;
    retp(y,os);
endp;

proc(2) = taildatan(x,l);
    local n,sdata,os,y;
    n=rows(x);
    sdata = sortc(x,1);
    os = sdata[l+1,1];
    y = os-sdata[1:l,1];
    retp(y,os);
endp;

/* xn are exceedances over osn */
l = round(n^(0.8-0.01));
{ xn,osn } = taildata(yn,l);
/*{ xn_n,osn_n } = taildatan(yn,l);*/

/*The next 11 lines solves the nonlinear equation two lines below equation (7) in Martins-Filho et al. (2018).*/
ystart = osn;
yleth = yn; 
afhatau = 1-l/n;
rothl = 0.79*iqr(yn)*n^(-1/5+0.01);
nlc.output = 0;
{nlc,nlo} = nlsys(nlc, &fhatau, ystart); 
osns = nlo.xp;
xns = yn - osns;
indixns = (xns .>zeros(n,1));
xns = selif(xns,indixns);
xns = sortc(xns,1);

/* Next is just a procedure to calculate interquartile range.*/
proc iqr(x);
    local s,i,n;
    n = rows(x);
    s = sortc(x,1);
    i = s[round(0.75*n)]-s[round(0.25*n)];
    retp(i);
endp;
/* Next is the function described in equation (7) in Martins-Filho et al. (2018).*/
proc fhatau(x);
    local k,sumf,fh,n,b,lb,onez,yleth,rothl;
    yleth = varget("yleth");
    rothl = varget("rothl");
    n = rows(yleth);
    k = (1/rothl)*(ones(n,1)*x-yleth);
    b = k;
    onez = ones(n,1);
    lb = -onez;
    sumf =  ( 0.75*(b-(1/3)*b^3 + (2/3)*onez) ).* ((-onez .<= b) .and (b .<= onez)) ;
    sumf = sumf + (b.>onez);
    fh = meanc(sumf) - afhatau;
    retp(fh);
endp;

/* L Moments Estimators of the GPD parameters.  Not used*/
parn = lmomgpd(xn);

/*print "parn" parn;*/
/*parn_n = lmomgpd(xn_n);*/

/* llgpd(par,x) is a pointer to the loglikelihood of a Generalized Pareto Distribution with scale and shape parameters. Location
parameter is taken to be zero. This will be used in a CML procedure. gradgpd(par,x) is a gradient procedure pointer for the CML.*/
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
        e = (0.< x);
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

/* Maximum Likelihood Estimation: use smooth quantile*/
/* Data structure */

struct DS d0;
d0 = reshape(dsCreate,1,1);
d0[1].dataMatrix = xns;
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
out1 = cmlmt(&llgpd,p1,d0,c1); 
bini = pvUnpack(out1.par,1);
phns = bini;
/*print "phns" phns;*/

proc(1) = zhath(par,n,os,q,l);
    local zhat;
    zhat = os*((n*(1-q)/l)^par);
    retp(zhat);
endp;

/* The next procedure calculates VaR and ES in accordance to equations (9) and (10) in Martins-Filho et al. (2018)*/
proc(2) = zhatq(par,n,os,q);
    local zhat,ez;
    zhat = os - (par[1]/par[2]).*( (n*(1-q)/l)^par[2] - 1 );
    ez = zhat*(1/(1+par[2])); 
    retp(zhat,ez);
endp;

/* Now we collect several quantiles (q): 0.95, 0.99, 0.995, 0.999*/

storeq = zeros(4,2);
alphvalue = {0.05,0.01,0.005,0.001};
alphi = 1; 
    do while alphi<=rows(alphvalue);
        alph = alphvalue[alphi];
        q=1-alph;   
        { mzhatns, mezns } = zhatq(phns,n,osns,q);
        /* The first entry is VaR and the second entry is ES*/
        storeq[alphi,.] = mzhatns~mezns; 
        alphi = alphi + 1;
    endo;
store[repeati,.] = phns'~osns~storeq[1,1]~storeq[2,1]~storeq[3,1]~storeq[4,1]~storeq[1,2]~storeq[2,2]~storeq[3,2]~storeq[4,2];
/*Elements of this row vector are:
 sigma (GPD) k(GPD) quantile(of order 1-l/n) 0.95-VaR 0.99-VaR 0.995-VaR 0.999-VaR 0.95-CES 0.99-CES 0.995-CES 0.999-CES
 The quantile, VaR and CES are for the error (residuals)*/

/* For quantile 0.95 */
storey1[repeati,.] = evy~evx~evmhat~sqrt(evshat2)~(evmhat+sqrt(evshat2)*storeq[1,1])'~(evmhat+sqrt(evshat2)*storeq[1,2])';
/* For quantile 0.99 */
storey2[repeati,.] = evy~evx~evmhat~sqrt(evshat2)~(evmhat+sqrt(evshat2)*storeq[2,1])'~(evmhat+sqrt(evshat2)*storeq[2,2])';
/* For quantile 0.995 */
storey3[repeati,.] = evy~evx~evmhat~sqrt(evshat2)~(evmhat+sqrt(evshat2)*storeq[3,1])'~(evmhat+sqrt(evshat2)*storeq[3,2])';
/* For quantile 0.999 */
storey4[repeati,.] = evy~evx~evmhat~sqrt(evshat2)~(evmhat+sqrt(evshat2)*storeq[4,1])'~(evmhat+sqrt(evshat2)*storeq[4,2])';
repeati = repeati + 1;
endo;
/* Print the results. */
/*store;
storey1;
storey2;
storey3;
storey4;*/

/* Backtest */

stsel = 1; 
edsel = repeatn;

/* Only for 95%, as only storey1 is used to define ytrue.*/
ytrue = storey1[stsel:edsel,1:4];

/* This is VaR for all quantile levels*/
xqn = storey1[stsel:edsel,5]~storey2[stsel:edsel,5]~storey3[stsel:edsel,5]~storey4[stsel:edsel,5];
/* This is ES for all quantile levels */
en = storey1[stsel:edsel,6]~storey2[stsel:edsel,6]~storey3[stsel:edsel,6]~storey4[stsel:edsel,6];

/* Bootstrap Test for ES */
indxqnd  = zeros(rows(alphvalue),1);
pvalxqnd = zeros(rows(alphvalue),1);
inxqnd   = zeros(rows(ytrue),rows(alphvalue));

prt = zeros(rows(alphvalue),1);

alphi = 1; 
do while alphi<=rows(alphvalue);
{indxqnd[alphi,1],pvalxqnd[alphi,1],inxqnd[.,alphi]}  = backt(ytrue[.,1],xqn[.,alphi],alphvalue[alphi]);

rt = (ytrue[.,1]- en[.,alphi])./ytrue[.,4];
rt = rt.*inxqnd[.,alphi];
rt = selif(rt,(rt ./= 0));
resampl=1000;
prt[alphi,1] = backtes(rt,resampl);
alphi = alphi + 1;
endo;

print "Results for our estimators";

print "1st: Backtest for n=" rows(ytrue) "Alpha=" alphvalue'; 
"The expected violation number is n*p= " (rows(ytrue)*alphvalue)';
print "2nd:actual violation number(exceedances over the quantile)";
indxqnd';
print "3rd is the p-value for VaR: compare to 5%";
pvalxqnd';
print "4th: The p-value for exceedance residuals in GPD have mean zero (CES):";
prt';

/*-----------------------------------------------------*/

@Backtest.g: x is the return series, and y is the qunatile, a=alph@
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

proc(1)= backtes(x,num);
    local n,tr,rbar,sigr,m,i,xstar,trs,rbars,sigrs,id,me,j;
    n=rows(x);
    rbar=(1/n)*sumc(x);
    sigr=sqrt((1/(n-1))*sumc((x-rbar)^2));
    tr=(rbar-0)/(sigr/(sqrt(n)));

    i=1;
    m=0;
    do while i<=num;
    id=n.*rndu(n,1);
    me=zeros(n,n);
        j=1;
        do while j<=n;
        me[.,j]=((j-1).<=id .and id .<j)*j;
        j=j+1;
        endo;
    id=sumc(me');
    id=id+n*((sumc(me')).==0);
    xstar=x[id];
    rbars=(1/n)*sumc(xstar);
    sigrs=sqrt((1/(n-1))*sumc((xstar-rbars)^2));
    trs=(rbars-rbar)/(sigrs/(sqrt(n)));
    if trs>tr;
    m=m+1;
    endif;
    i=i+1;
    endo;
    retp(m/num);
endp;

struct plotcontrol myplot;
myplot = plotGetDefaults("xy");
t = seqa(1,1, rows(ytrue));
plotSetXLabel(&myplot,"time");
plotSetYLabel(&myplot,"log returns");

_pltype = { 6 };
gr = sortc(t~storey1[.,1]~storey1[.,5]~storey1[.,6],1);
plotxy(myplot,gr[.,1],gr[.,2]~gr[.,3]~gr[.,4]);
plotSave("D:\\projects\\IFPRI\\FSP\\gauss\\outputfolder\\corn_beta.jpg", 30|18);
/*plotSave()*/

n1=rows(gr);
period=60; 
corvec=zeros(n1-period+1,1); 
indvec=zeros(n1-period+1,1); 
pvec=zeros(n1-period+1,1);

l=1;
do while l<=n1-period+1;
{indr,pvalr,inr,corx}=backtc(gr[.,2],gr[.,2],period,l,0.95);
corvec[l]=corx;
indvec[l]=indr;
pvec[l]=pvalr;
l=l+1;
endo;

alpha=0.05;
vecviol = ones(n1-period+1,1)*period*alpha;
print "First column is a time index, second column are returns, third column is the 95 percent conditional quantile fourth column is the color code";
print gr[period:n1,1]~gr[period:n1,2]~gr[period:n1,3]~gr[period:n1,4]~corvec;/*change*/



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


















proc mupx(y,x,z,p,h,h0,pt,vpx);
    local zeval,n,b,a,k,by,ay,aux,f,mup;
    zeval = z;
    n = rows(y);
    b=x-zeval*ones(rows(x),1);
    a=inv(h)*b;
    k=(abs(a).<=ones(n,1)).*(0.75*(ones(n,1)-a^2));
    by = vpx*ones(n,1)-y;
    ay = inv(h0)*by;
    aux = y.*(ones(n,1) -cdfepa(ay)) +h0*evepa(ay);
    f = inv(n*h)*sumc(k.*pt);
    mup = inv(p*n*h)*sumc((k.*pt).*aux)./f;
    retp(mup);
endp;

proc fc(x);
    local y,z,zeval,h0,h,p,n,pt,b,a,k,by,ay,aux,f,cdf_wdkll;
    y = varget("y");
    z = varget("x");
    zeval = varget("xeval");
    h0 = varget("hy");
    h = varget("hx");
    p = varget("peval");
    n = rows(y);
    pt = varget("pteval");
    b=z-zeval*ones(rows(z),1);
    a=inv(h)*b;
    k=(abs(a).<=ones(n,1)).*(0.75*(ones(n,1)-a^2));
    by = (x*ones(n,1)-y);
    ay = inv(h0)*by;
    aux = cdfepa(ay);
    f = inv(n*h)*sumc(k.*pt);
    cdf_wdkll = inv(n*h)*sumc((k.*pt).*aux)./f;
    retp(cdf_wdkll-p);
endp;

proc cdfepa(x);
	local g,n,m;
    n = rows(x);
    m = cols(x);
    g = (x.>-1*ones(n,m)).*(0.75*(-x^3/3+x+2/3));
    g = g.*(x.<=ones(n,m))+(x.>ones(n,m));
    retp(g);
endp;

proc evepa(x);
    local g,n,m;
    n = rows(x);
    m = cols(x);
    g = (x.>=-1*ones(n,m)).*(0.75*(0.25-0.5*x^2+0.25*x^4));
    g = g.*(x.<=ones(n,m));
    retp(g);
endp;

proc aicce(x);
	local result,y,z,vp,ys,n,pt,b,a,k,wi,hc,sig2,th,h0,by,ay;
    y = varget("y");
    z = varget("x");
    vp = varget("xqnca");
    h0 = varget("hy");
    n = rows(y);
    by = (vp*ones(n,1)-y);
    ay = inv(h0)*by;
    ys = y.*(ones(n,1) -cdfepa(ay)) +h0*evepa(ay);
    /* pt = caivpt(y,z,x); */
    pt = varget("ptt");
    b=z-z';
    a=inv(x)*b;
    k=(abs(a).<=ones(n,n)).*(0.75*(ones(n,n)-a^2));
    k = inv(x)*k;
    wi = ((sumc(pt.*k))^(-1))*ones(1,n);
    hc = (pt.*k)'.*wi;
    sig2 = (eye(n) - hc)*ys;
    sig2 = (1/n)*moment(sig2,0);
    th = sumc(diag(hc));
    result = ln(sig2)+(n+th)/(n-(th+2));
    retp(result);
endp;

proc aiccv(x);
	local result,y,z,vp,ys,n,pt,b,a,k,wi,hc,sig2,th,h0,by,ay;
    y = varget("y");
    z = varget("x");
    vp = varget("xqnca");
    h0 = varget("h0");
    n = rows(y);
    by = (vp*ones(n,1)-y);
    ay = inv(h0)*by;
    ys = cdfepa(ay);
    /* pt = caivpt(y,z,x); */
    pt = varget("ptt");
    b=z-z';
    a=inv(x)*b;
    k=(abs(a).<=ones(n,n)).*(0.75*(ones(n,n)-a^2));
    k = inv(x)*k;
    wi = ((sumc(pt.*k))^(-1))*ones(1,n);
    hc = (pt.*k)'.*wi;
    sig2 = (eye(n) - hc)*ys;
    sig2 = (1/n)*moment(sig2,0);
    th = sumc(diag(hc));
    result = ln(sig2)+(n+th)/(n-(th+2));
    retp(result);
endp;


proc aiccf(x);
	local result,y,z,vp,ys,n,pt,b,a,k,wi,hc,sig2,th;
    y = varget("y");
    z = varget("x");
    vp = varget("xqnca");
    n = rows(y);
    ys = (y .<= ones(n,1)*vp);
    /* pt = caivpt(y,z,x); */
    pt = varget("ptt");
    b=z-z';
    a=inv(x)*b;
    k=(abs(a).<=ones(n,n)).*(0.75*(ones(n,n)-a^2));
    k = inv(x)*k;
    wi = ((sumc(pt.*k))^(-1))*ones(1,n);
    hc = (pt.*k)'.*wi;
    sig2 = (eye(n) - hc)*ys;
    sig2 = (1/n)*moment(sig2,0);
    th = sumc(diag(hc));
    result = ln(sig2)+(n+th)/(n-(th+2));
    retp(result);
endp;

proc(1) = ullg(y,x,h,z);
    local n,m,b,a,k,a1,a2,a3,b1,b2;
    n = rows(y); 
    m = zeros(rows(z),1); 
    b = (x-z');
    a = inv(h)*b;
    k = (abs(a).<=ones(n,cols(b))).*(0.75*(ones(n,cols(b))-a^2));
    a1 = sumc(k); a2 = sumc(k.*b); a3 = sumc(k.*b^2);
    b1 = sumc(k.*y); b2 = sumc((k.*b).*y);
    m = (a3.*b1 - a2.*b2)./(a1.*a3-a2^2+ n^(-2));
    retp(m);
endp;

proc(1)=mllgx(y,x,h,z);

    local n,e,d,m_ll,mh,b,a,k,w,aux,j,m_lls,e1,dm;

    n=rows(y);d=cols(x);m_ll=zeros(rows(z),1); e=1~zeros(d,1)';
    mh=diagrv(zeros(d,d),h'); e1 = (zeros(d,1)')~1; m_lls = zeros(rows(z),1);

    if d==1;
        dm = 0|n^(-2);
    elseif d==2;
        dm = 0|n^(-2)|n^(-2);
    endif;
    dm = diagrv(zeros(d+1,d+1),dm);

    j=1; 
    do while j<=rows(z);
    b=(x'-ones(n,1)'.*.z[j,.]');
    a=inv(mh)*b;
    @k = (abs(a).<=sqrt(5)).*((0.75/sqrt(5))*(ones(d,n)-0.2*a^2));@
    k=(abs(a).<=1).*(0.75*(ones(d,n)-a^2));
    @k=inv(sqrt(2*pi))*exp(-0.5*a^2);@
    w=diagrv(zeros(n,n),prodc(k));
    aux=ones(n,1)~b';
    m_ll[j]=e*inv(aux'*w*aux+dm)*aux'*w*y;
    @m_lls[j] = e1*inv(aux'*w*aux)*aux'*w*y;@
    j=j+1;
    endo;
    retp(m_ll);
endp;

proc(1) = mllg(y,x,h,z);

    local n,e,d,m_ll,mh,j,b,a,k,w,aux;

    n = rows(y); d = cols(x);m_ll = zeros(rows(z),1); e = 1~zeros(d,1)';
    mh = diagrv(zeros(d,d),h); 

    j = 1;
    do while j <= rows(z);
    b = (x'-ones(n,1)'.*.z[j,.]');
    a = inv(mh)*b;
    k = (abs(a).<=1).*(0.75*(ones(d,n)-a^2));
    w = diagrv(zeros(n,n),prodc(k));
    aux = ones(n,1)~b';
    m_ll[j] = e*inv(aux'*w*aux)*aux'*w*y;
    j = j+1;
    endo;
    retp(m_ll);

endp;

proc(1)= mnwgx(y,x,h,z);

    local n,d,m_nw,mh,b,a,k,w,aux,j,m_nwf,f;

    n=rows(y);d=cols(x);m_nw=zeros(rows(z),1); 
    m_nwf = zeros(rows(z),1); f = zeros(rows(z),1);
    mh=diagrv(zeros(d,d),h); 

    j=1; 
    do while j<=rows(z);
    b=(x'-ones(n,1)'.*.z[j,.]');
    a=inv(mh)*b;
    @k = (abs(a).<=1).*((15/32)*(7*a^4 -10*a^2-3*(ones(d,n))));@
    k=(abs(a).<=1).*(0.75*(ones(d,n)-a^2));
    @k=inv(sqrt(2*pi))*exp(-0.5*a^2);@
    w=diagrv(zeros(n,n),prodc(k));
    aux=ones(n,1);
    f[j] = inv(n*prodc(h))*aux'*w*aux;
    m_nwf[j] = inv(n*prodc(h))*aux'*w*y; 
    if f[j].>0;
    m_nw[j]=inv(aux'*w*aux)*aux'*w*y;
    elseif f[j].==0;
    m_nw[j]=m_nwf[j]/(f[j]+1/n);	
    endif;
    j=j+1;
    endo;
    retp(m_nw);
endp;

proc(1)= nw(y,x,h);
    local n,m_nw,b,a,k,aux,f;
    n=rows(y);m_nw=zeros(n,1); 
    b=(x-x');
    a=inv(h)*b;
    k=(abs(a).<=ones(n,n)).*(0.75*(ones(n,n)-a^2));
    @k=inv(sqrt(2*pi))*exp(-0.5*a^2);@
    f = inv(n*h)*sumc(k);
    aux = y*ones(1,n);
    m_nw = inv(n*h)*sumc(k.*aux)./f;
    retp(m_nw);
endp;

proc(1)= nwe(y,x,h,z);
    local n,m_nw,b,a,k,aux,f,m;
    n=rows(y);m_nw=zeros(n,1);
    m = rows(z);
    b=(x-z');
    a=inv(h)*b;
    k=(abs(a).<=ones(n,m)).*(0.75*(ones(n,m)-a^2));
    @k=inv(sqrt(2*pi))*exp(-0.5*a^2);@
    f = inv(n*h)*sumc(k);
    aux = y*ones(1,m);
    if f.==zeros(m,1);
    m_nw = inv(n*h)*sumc(k.*aux)./(f+1/n);
    else;
    m_nw = inv(n*h)*sumc(k.*aux)./f;
    endif;
    retp(m_nw);
endp;

proc(1) = f_hat(x,h,z);
    local a,k,f,n,b,m;
    n = rows(x); 
    m = rows(z);
    b = x - z';
    a = inv(h)*b;
    k = (abs(a).<=ones(n,m)).*(0.75*(ones(n,m)-a^2));
    f = inv(n*h)*sumc(k);
    retp(f);
endp;

/********************************************************************
The procedure LAMBDA_OPT has the format 
{lambda,f,p,ret}=lambda_opt(e,output);

Input:
e		nxm 	moment conditions
output		1x1	1 if iteration output printed to screen

Output:
lambda		mx1	final value
f		1x1	log-likelihood
p		nx1	probability vector
ret		1x1	return code
			ret=0 for normal convergence
			ret=1 if maximum iterations reached
			ret=2 if failure to find a stepsize
********************************************************************/

proc (4) = lambda_opt(e,_output);
local n,n1,maxiters,maxstepiters,maxtol,lambda,w,m,em,
lik,u,direc,step,lambda1,lik1,p,ret,iv,ev,eh,test,uu,su,
direc0,it,si,lambda2,lik2,steps,mineig,lambda3,lik3;  
  MaxIters=100; 		@ maximum number of newton iterations 	@
  MaxStepIters=13;		@ maximum number of steplength steps  	@
  MaxTol = .00001;		@ tolerance for convergence		@
  MinEig =.1;			@ minimum eigenvalue of Hessian		@
  n=rows(e);
  n1=1/n;
  lambda=zeros(cols(e),1);
  w=ones(n,1);
  ret=0;
  lik=0;
  uu=moment(e,0);
  trap 1;
    m=inv(chol(uu));
  trap 0;
  if scalerr(m);
    {ev,eh}=eighv(uu);
    iv=(ev.*(ev .> mineig) + mineig .* (ev .<= mineig));
    m=eh.*(1./sqrt(iv));
  endif;
  em=e*m;
  for it (1,maxiters,1);
    u=em./w;
    uu=moment(u,0);
    su=sumc(u);
    trap 1;
      direc=su/uu;
    trap 0;
    if scalerr(direc);
      {ev,eh}=eighv(uu);
      iv=1./(ev.*(ev .> mineig) + mineig .* (ev .<= mineig));
      direc=(eh*(iv.*(eh')))*su;
    endif;
    step=1;  
    for si (1,maxstepiters,1);
      steps=si;
      lambda1=lambda+direc*step;  
      w=1+em*lambda1; 
      if (minc(w) .> n1);
        lik1=-sumc(ln(w));
        if lik1 < lik; 
          lambda2=lambda+direc*step/2;
          lik2=-sumc(ln(1+em*lambda2));
          break;
        endif;
      endif;
      step=step/2;
    endfor;
    if steps==maxstepiters; ret=2; break; endif;
    step=(lik1+3*lik-4*lik2)/(4*lik1+4*lik-8*lik2);
    lambda3=lambda+direc*step;  
    if (minc(1+em*lambda3) .> n1); 
      lik3=-sumc(ln(1+em*lambda3));
      si=minindc(lik1|lik2|lik3);
      lambda=lambda1*(si==1)+lambda2*(si==2)+lambda3*(si==3);
    else;
      si=minindc(lik1|lik2);
      lambda=lambda1*(si==1)+lambda2*(si==2);
    endif; 
    w=1+em*lambda; 
    lik=-sumc(ln(w));
    test=maxc(abs(em'(1./w)))/n;
    if _output==1;
      format 16,8;
      "=======================";
      "iteration  " it;
      "Steps      " steps;
      "Steplength " step;
      "function   " lik;
      "Max Error  " test;
      "_______________________";
      "  Parameter        Value ";
      seqa(1,1,rows(lambda))~(m*lambda);
      "";"";
    endif;
    if steps==maxstepiters; ret=2; break; endif;
    if (test < MaxTol); break; endif; 
    if it==maxiters; ret=1; endif;
  endfor;
  p=1./((1+em*lambda)*n);
  p=p./sumc(p);
retp(m*lambda,-lik,p,ret);
endp;

