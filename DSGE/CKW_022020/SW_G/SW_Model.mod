// Smets-Wouters Model

var ewma epinfma yf y r a b g qs ms spinf sw kpf kp cf invef c inve pinf w lab zcapf rkf kf pkf labf wf rrf mc zcap rk k pk;    
 
varexo eg em ea eb eqs epinf ew;  
 
parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
    czcap csadjcost ctou csigma chabb cfc 
    cindw cprobw cindp cprobp csigl clandaw 
    crpi crdy cry crr 
    crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow  
    ctrend cg
    sigma_a sigma_b sigma_g sigma_qs sigma_m sigma_pinf sigma_w;

// fixed parameters
ctou=.025;
clandaw=1.5;
cg=0.18;
curvp=10;
curvw=10;

// estimated parameters, mode
calfa	     = 0.2024;
cbeta	     = 0.9995;
csigma	     = 1.2679;
cfc	     = 1.6670;
cgy	     = 0.5881;
csadjcost    = 6.3144;
chabb	     = 0.8056;    
cprobw	     = 0.7668;
csigl	     = 2.5201;
cprobp	     = 0.5304;
cindw	     = 0.5345;
cindp	     = 0.1779;
czcap	     = 0.3597;
crpi	     = 1.8685;
crr	     = 0.8739;
cry	     = 0.1203;
crdy	     = 0.1282;
crhoa	     = 0.9826;
crhob	     = 0.1391;
//crhog	     = 0.9686;
crhog	     = 0.7;
crhols	     = 0.9928;
crhoqs	     = 0.6121;
crhoas	     = 1;
crhoms	     = 0.1999;
crhopinf     = 0.9856;
crhow	     = 0.9818;
cmap 	     = 0.8340;
cmaw  	     = 0.9337;
constelab    = 1.3263;
constepinf   = 0.6365;
constebeta   = 0.1126;
ctrend       = 0.5113;

sigma_a = 0.5017;
sigma_b = 0.3583;
sigma_g = 0.6752;
sigma_qs = 0.5678;
sigma_m = 0.2290;
sigma_pinf = 0.2181;
sigma_w = 0.2290;

// model solution

model(linear);

// parameter relations

#cpie=1+constepinf/100;
#cgamma=1+ctrend/100 ;
#cbeta=1/(1+constebeta/100);

#clandap=cfc;
#cbetabar=cbeta*cgamma^(-csigma);
#cr=cpie/(cbeta*cgamma^(-csigma));
#crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
#cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
//cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*((cbeta^(-1))*(cgamma^csigma) - (1-ctou))^calfa))^(1/(1-calfa));
#cikbar=(1-(1-ctou)/cgamma);
#cik=(1-(1-ctou)/cgamma)*cgamma;
#clk=((1-calfa)/calfa)*(crk/cw);
#cky=cfc*(clk)^(calfa-1);
#ciy=cik*cky;
#ccy=1-cg-cik*cky;
#crkky=crk*cky;
#cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
#cwly=1-crk*cky;

#conster=(cr-1)*100;

// flex-price benchmark economy

0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf);

zcapf =  (1/(czcap/(1-czcap)))* rkf;

rkf =  (wf)+labf-kf;

kf =  kpf(-1)+zcapf;

invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf )+qs;

pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1);

cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b;

yf = ccy*cf+ciy*invef+g  +  crkky*zcapf;

yf = cfc*( calfa*kf+(1-calfa)*labf +a);

wf = csigl*labf +(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1);

kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs;

// sticky price/wage economy

mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a;

zcap =  (1/(czcap/(1-czcap)))* rk;

rk =  w+lab-k;

k =  kp(-1)+zcap;

inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs;

pk = -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1);

c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b;

y = ccy*c+ciy*inve+g  +  1*crkky*zcap;

y = cfc*( calfa*k+(1-calfa)*lab +a );

pinf = (1/(1+cbetabar*cgamma*cindp))*(cbetabar*cgamma*pinf(1)+cindp*pinf(-1)               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc))+spinf;
 
w =  (1/(1+cbetabar*cgamma))*w(-1) +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
+(cindw/(1+cbetabar*cgamma))*pinf(-1)-(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
+(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)+(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*(csigl*lab + (1/(1-chabb/cgamma))*c -((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w)+ 1*sw;

r =  crpi*(1-crr)*pinf+cry*(1-crr)*(y-yf)+crdy*(y-yf-y(-1)+yf(-1))+crr*r(-1)+ms;

a = crhoa*a(-1)  + sigma_a * ea;

b = crhob*b(-1) + sigma_b * eb;

g = crhog*(g(-1)) + sigma_g * eg + cgy*sigma_a*ea;

qs = crhoqs*qs(-1) + sigma_qs * eqs;

ms = crhoms*ms(-1) + sigma_m * em;

spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);

epinfma = sigma_pinf * epinf;

sw = crhow*sw(-1) + ewma - cmaw*ewma(-1);

ewma = sigma_w * ew;
 
kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs;

end; 

shocks;
var ea;
stderr 1;
var eb;
stderr 1;
var eg;
stderr 1;
var eqs;
stderr 1;
var em;
stderr 1;
var epinf;
stderr 1;
var ew;
stderr 1;
end;

stoch_simul(order=1,periods=200,irf=0);