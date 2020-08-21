//Pertenece a Rodson Jeffrey Bejarano Torres
// Docente: Prof. Dr. Juan C. Cutipa Luque 
//Lista de Ejercicios 3

clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows

// Parametros de la planta /////
kr = 0.02
kb = 0.02
r1 = 0.015;
c2 = kr/r1;
c1 = kr*kb/(r1)^2;
g = 9.81;
m1 = 0.1;
m2 = 1;
l1=0.25;
// Modelo de la planta //
ap = [0.001 0 1 0; 0 0 0 1; 0 -m1*g/m2 -c1/m2 0; 0 (m1+m2)*g/(m2*l1) c1/(m2*l1) 0];
bp = [0 0; 0 0; c2/m2 0;  0 -c2/(m2*l1) ];
cp = [1 0 0 0; 0 1 0 0];
dp = 0*ones (2,2);

// Controlabilidad y Observabilidad
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap,bp)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap, cp)
rankO=rank(O)

ap_=ap
bp_=bp
cp_=cp
dp_=dp

// Se aumenta a la planta integradores antes del proyecto de controlador

[ns,nc] = size(bp_);                     //ns = número de entradas;  
                                        //nc = número de controles;   
a_1 = [ap_            bp_ ; 
       0*ones(nc,ns)  0*ones(nc,nc) ];
b_1 = [0*ones(ns,nc); eye(nc,nc)];
c_1 = [cp_      0*ones(nc,nc)];
d_1 = 0*ones(nc,nc)

G= syslin("c",a_1, b_1, c_1, d_1)
w = logspace(-3,3);

ms=1.7;// 0.3;%1.5;    % guarantee overshot Mp < 6dB = 20*log10(2) 
wbs=0.23;//0.05;%0.23;
ee=1e-3;//1e-4
ki=1; // used to give more accurate adjustment to the cut-off frequency wbs
      // by default set it to 1
//           --------     WT Data    ------------
mt=1.3;//1.00;    % guarantee overshot Mp < 2dB = 20*log10(1.26)
wbt=4.1;//9.1;%4.1;
ee=1e-3;//1e-4


//  WS   

s=poly(0,'s');
ws1=(s/ms+wbs)/(s+wbs*ee),
ws2=ws1;
ws=[ws1,0;0,ws2]
//Ws=syslin('c',ws)
Ws=blockdiag(ws1,ws2)

//   WT   

s=poly(0,'s');
wt1=(s+wbt/mt)/(ee*s+wbt),
wt2=wt1;
wt=[wt1,0;0,wt2]
//Wt=syslin('c',wt)
Wt=blockdiag(wt1,wt2)


//   WR   
s=poly(0,'s');
wr1=s/s,
wr2=wr1;
wr=[wr1,0;0,wr2]



// ------------------ Funciones de ponderación de trazado

svs = svplot(Ws,w);
svt = svplot(Wt,w);
scf(2);
plot2d("ln", w, [-20*log(svs')/log(10) -20*log(svt')/log(10)])
xgrid(12)
xtitle("Singular values plot inv(Ws) and inv(Wt)","Frequency (rad/s)", "Amplitude (dB)");


[P,r]=augment(G,'ST');
//[P,r]=augment(g,'SRT');
P = blockdiag(Ws,Wt,eye(G))*P;
//P=minreal(P);

//trick to tackle when "D12 is not full rank"
P.D(1,3)=0.0001;
P.D(2,4)=0.0001;

r=[2,2]
romin=0.0001
romax=2000
nmax=100;
//[K,ro]=h_inf(P,r,romin,romax,nmax)
// alternatives
//[AK,BK,CK,DK,(RCOND)] = hinf(P.A,P.B,P.C,P.D,2,2,4)
K = ccontrg(P, [2,2], 1.3) // this is good for me and for this system

// -------------- Analisis de el  Feeedback Control System

[Se,Re,Te]=sensi(G,K) // S=(I+GK)^-1, T=I-S=GK(I+GK)^-1

// ------------------ Plot weighting functions

svS = svplot(Se,w);
svT = svplot(Te,w);
scf(3);
plot2d("ln", w,[-20*log(svs')/log(10) 20*log(svS')/log(10)],[-1 -1 2 2],leg="$\overline{\sigma}(W_S^{-1})$@$\underline{\sigma}(W_S^{-1})$@$\overline{\sigma}(S)$@$\underline{\sigma}(S)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
//set(gca(),"auto_clear","off")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");


scf(4);
plot2d("ln", w,[-20*log(svt')/log(10) 20*log(svT')/log(10)],[-1 -1 2 2],leg="$\overline{\sigma}(W_T^{-1})$@$\underline{\sigma}(W_T^{-1})$@$\overline{\sigma}(T)$@$\underline{\sigma}(T)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
xtitle("","Frequency (rad/s)", "Amplitude (dB)");

// --------------- Open loop and Closed loop analysis --------------

sysOL=G*K
sysCL=G*K*inv(eye(2,2)+G*K)
// eigenvalues in LHP (Left Half Plane), stable according to RH
spec(sysCL.A)

sv1= svplot(sysOL,w);
sv2= svplot(sysCL,w);
scf(5);
plot2d("ln", w,[20*log(sv1')/log(10) 20*log(sv2')/log(10)],[2 2 3 3],leg="$\overline{\sigma}(GK)$@$\underline{\sigma}(GK)$@$\overline{\sigma}(T)$@$\underline{\sigma}(T)$")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");
xgrid(12)
//set(gca(),"auto_clear","off")
xtitle("","Frequency (rad/s)", "Amplitude (dB)");

// Time responses in XCOS
