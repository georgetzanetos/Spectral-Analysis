 clear all
clc 

%%
% Desig and Turbulence Parameters 
V   = 51.4;
Spw   = 24.2;
b   = 13.36;
mub = 13;
KX2 = 0.016;
KZ2 = 0.040;
KXZ = 0.002;
CL  = 1.1360;

W         = 59143;
rho       = 1.225;

Lg        = 150; 
B         = b/(2*Lg);
sigma     = 1;
sigmaug_V = sigma/V;
sigmavg   = sigma;
sigmabg   = sigmavg/V;
sigmaag   = sigma/V;

Iug0 = 0.0249*sigmaug_V^2;
Iag0 = 0.0182*sigmaag^2;
tau1 = 0.0991;     tau2 = 0.5545;     tau3 = 0.4159;
tau4 = 0.0600;     tau5 = 0.3294;     tau6 = 0.2243;


CYb  =-0.9896;     Clb  =-0.0772;     Cnb  = 0.1628;
CYp  =-0.0870;     Clp  =-0.3415;     Cnp  =-0.0268;
CYr  = 0.4300;     Clr  = 0.2830;     Cnr  =-0.1930;
CYda = 0.0000;     Clda =-0.2349;     Cnda = 0.0286;
CYdr = 0.3037;     Cldr = 0.0286;     Cndr =-0.1261;
 
                   Clpw = 0.8*Clp;    Cnpw = 0.9*Cnp;
                   Clrw = 0.7*Clr;    Cnrw = 0.2*Cnr;
CYfb = 0;
Clfb = 0;
Cnfb = 0;


%%
% Matrix Entries
yb   = (V/b)*CYb/(2*mub);
yphi = (V/b)*CL/(2*mub);
yp   = (V/b)*CYp/(2*mub);
yr   = (V/b)*(CYr-4*mub)/(2*mub);
ybg  = yb;
ydr  = (V/b)*CYdr/(2*mub);
den  = b*4*mub*(KX2*KZ2-KXZ^2)/V;
lb   = (Clb*KZ2+Cnb*KXZ)/den;
lp   = (Clp*KZ2+Cnp*KXZ)/den;
lr   = (Clr*KZ2+Cnr*KXZ)/den;
lda  = (Clda*KZ2+Cnda*KXZ)/den;
ldr  = (Cldr*KZ2+Cndr*KXZ)/den;
lug  = (-Clrw*KZ2-Cnrw*KXZ)/den;
lbg  = lb;
lag  = (Clpw*KZ2+Cnpw*KXZ)/den;
nb   = (Clb*KXZ+Cnb*KX2)/den;
np   = (Clp*KXZ+Cnp*KX2)/den;
nr   = (Clr*KXZ+Cnr*KX2)/den;
nda  = (Clda*KXZ+Cnda*KX2)/den;
ndr  = (Cldr*KXZ+Cndr*KX2)/den;
nug  = (-Clrw*KXZ-Cnrw*KX2)/den;
nbg  = nb;
nag  = (Clpw*KXZ+Cnpw*KX2)/den;
aug1 =-(V/Lg)^2*(1/(tau1*tau2));
aug2 =-(tau1+tau2)*(V/Lg)/(tau1*tau2);
aag1 =-(V/Lg)^2*(1/(tau4*tau5));
aag2 =-(tau4+tau5)*(V/Lg)/(tau4*tau5);
abg1 =-(V/Lg)^2;
abg2 =-2*(V/Lg);
bug1 = tau3*sqrt(Iug0*V/Lg)/(tau1*tau2);
bug2 = (1-tau3*(tau1+tau2)/(tau1*tau2))*sqrt(Iug0*(V/Lg)^3)/(tau1*tau2);
bag1 = tau6*sqrt(Iag0*V/Lg)/(tau4*tau5);
bag2 = (1-tau6*(tau4+tau5)/(tau4*tau5))*sqrt(Iag0*(V/Lg)^3)/(tau4*tau5);
bbg1 = sigmabg*sqrt(3*V/Lg);
bbg2 = (1-2*sqrt(3))*sigmabg*sqrt((V/Lg)^3);

%%
% Uncontrolled System
A = [yb yphi yp    yr 0    0    0    0    ybg  0;
     0  0    2*V/b 0  0    0    0    0    0    0;
     lb 0    lp    lr lug  0    lag  0    lbg  0;
     nb 0    np    nr nug  0    nag  0    nbg  0;
     0  0    0     0  0    1    0    0    0    0;
     0  0    0     0  aug1 aug2 0    0    0    0;
     0  0    0     0  0    0    0    1    0    0;
     0  0    0     0  0    0    aag1 aag2 0    0;
     0  0    0     0  0    0    0    0    0    1;
     0  0    0     0  0    0    0    0    abg1 abg2];

B = [0   ydr 0    0    0;
     0   0   0    0    0;
     lda ldr 0    0    0;
     nda ndr 0    0    0;
     0   0   bug1 0    0;
     0   0   bug2 0    0;
     0   0   0    bag1 0;
     0   0   0    bag2 0;
     0   0   0    0    bbg1;
     0   0   0    0    bbg2];
 
C = [180/pi 0      0     0     0 0 0 0 0 0 ;
     0      180/pi 0     0     0 0 0 0 0 0 ;
     0      0      2*V/b 0     0 0 0 0 0 0 ;
     0      0      0     2*V/b 0 0 0 0 0 0];

D = [0 0 0 0 0 ;
     0 0 0 0 0 ;
     0 0 0 0 0 ;
     0 0 0 0 0];
 
sys = ss(A,B,C,D);

pole(sys)


pzmap(sys)
title('Pole-Zero Map of Uncontrolled System')  
hm = findobj(gca, 'Type', 'Line');                        
hm(3).MarkerSize = 12; 
hm(3).LineWidth = 1.5; 
grid on
pause

% Checking for controllability 
ctrb_rank = rank(ctrb(sys.A,sys.B)); %full rank so controllable

% Feedback gain 
Kphi = -0.1;
K = [0 Kphi 0 0 0 0 0 0 0 0];
Acon = A-B(:,1)*K;
newsys = ss(Acon,B,C,D);

pole(newsys)

pzmap(newsys)
title('Pole-Zero Map of Controlled System')
gm = findobj(gca, 'Type', 'Line');                        
gm(3).MarkerSize = 12; 
gm(3).LineWidth = 1.5; 
grid on
pause

%% 
P = [(1/2)*b/V      0           0                   0;
    0           -(1/2)*b/V      0                   0;
    0               0      -4*mub*KX2*b/V       4*mub*KXZ*b/V;
    0               0      4*mub*KXZ*b/V       -4*mub*KZ2*b/V];

Q = [0              0           0                   -1;
     0              0          -1                  0;
    -Clb             0          -Clp                -Clr;
    -Cnb             0          -Cnp                -Cnr];

R = [0      0; 
     0      0;
    -Clda -Cldr;
    -Cnda -Cndr];


Ashort = inv(P)*Q;
Bshort = inv(P)*R;


Ashort2 = [0    0    0    0    0    0; 
           0    0    0    0    0    0;
          lug   0    lag  0    lbg  0;
          nug   0    nag  0    nbg  0;];
      
Ashort3 = [0  0    0     0  0    1    0    0    0    0;
           0  0    0     0  aug1 aug2 0    0    0    0;
           0  0    0     0  0    0    0    1    0    0;
           0  0    0     0  0    0    aag1 aag2 0    0;
           0  0    0     0  0    0    0    0    0    1;
           0  0    0     0  0    0    0    0    abg1 abg2];

Afull = [Ashort   Ashort2;
         Ashort3        ];
    

Bshort2 = [0    0    0;
           0    0    0;
           0    0    0;
           0    0    0];  
     
    
Bshort3 = [0   0   bug1  0    0;
           0   0   bug2  0    0;
           0   0   0    bag1 0;
           0   0   0    bag2 0;
           0   0   0    0    bbg1;
           0   0   0    0    bbg2];

Bred = [Bshort Bshort2;
        Bshort3      ];  
       
       
Cred = [180/pi    0       0         0         0   0  0  0    0    0;
        0      180/pi     0         0         0   0  0  0    0    0;
        0       0     2*V/b         0         0   0  0  0    0    0;  
        0       0     0         2*V/b         0   0  0  0    0    0;
        V*yb    V*yphi   V*yp  V*yr+V*2*V/b   0   0  0  0   V*ybg  0]; % a_y
     

 
Dred = [0  0    0     0  0;
        0  0    0     0  0;
        0  0    0     0  0;
        0  0    0     0  0;
        0  ydr  0    0   0];% a_y

  
% Feedback gain 2
Kphi2 = -0.1;
K2    = [0 Kphi2 0 0  0 0  0 0  0 0];
Ared   = Afull-Bred(:,1)*K2;
redsys = ss(Ared,Bred,Cred,Dred);

pole(redsys)

pzmap(redsys)
title('Pole-Zero Map of Reduced System')
gm = findobj(gca, 'Type', 'Line');                        
gm(3).MarkerSize = 12; 
gm(3).LineWidth = 1.5; 
figure(1)
grid on
pause

%% Time Domain Simulation
dt = 0.05; T  = 60; t = [0:dt:T]; N = length(t);
nn = zeros(1,N);

% turbulence inputs
v_g = randn(1,N)/sqrt(dt);
w_g = randn(1,N)/sqrt(dt);

w1 = [nn' nn' nn' nn'  nn'];     
w2 = [nn' nn' nn'  w_g' nn'];
w3 = [nn' nn' nn'  nn'  v_g'];

%% 
% Controlled System
Cexpanded = [180/pi    0       0         0         0   0  0  0    0    0;
             0      180/pi     0         0         0   0  0  0    0    0;
             0       0     2*V/b         0         0   0  0  0    0    0;  
             0       0     0         2*V/b         0   0  0  0    0    0;
             V*yb    V*yphi   V*yp  V*yr+V*2*V/b   0   0  0  0  V*ybg  0]; % a_y
         
Dexpanded = [0 0   0 0 0 ;
             0 0   0 0 0 ;
             0 0   0 0 0 ;
             0 0   0 0 0 ;
             0 ydr 0 0 0]; %a_y

% response to u_g
y1 = lsim(Acon,B,Cexpanded,Dexpanded,w1,t);
% response to v_g
y2 = lsim(Acon,B,Cexpanded,Dexpanded,w2,t);
% response to w_g
y3 = lsim(Acon,B,Cexpanded,Dexpanded,w3,t);
% response to combined set (superposition)
yt = y1+y2+y3;

%%
% Reduced System

% response to u_g
yy1 = lsim(Ared,Bred,Cred,Dred,w1,t);
% response to v_g
yy2 = lsim(Ared,Bred,Cred,Dred,w2,t);
% response to w_g
yy3 = lsim(Ared,Bred,Cred,Dred,w3,t);
% response to combined set (superposition)
yyt = yy1+yy2+yy3;

%% Plotting


subplot(5,1,1)
plot(t,[yt(:,1) yyt(:,1)]); 
% axis([0 60 -0.1  0.1]); 
xlabel('time [s]'); ylabel('beta [deg]');
pbaspect([5 1 1]);
lgd = legend('Controlled','Reduced')
legend('Location','north')
lgd.FontSize = 7

% figure(2)
subplot(5,1,2)
plot(t,[yt(:,2) yyt(:,2)]);
% axis([0 60 -0.2  0.2]); 
xlabel('time [s]'); ylabel('phi [deg]');
pbaspect([5 1 1]);

% figure(3)
subplot(5,1,3)
plot(t,[yt(:,3) yyt(:,3)]); 
% axis([0 60 -1.5e-2  1.5e-2]); 
xlabel('time [s]'); ylabel('pb/2V [deg]');
pbaspect([5 1 1]);

% figure(4)
subplot(5,1,4)
plot(t,[yt(:,4) yyt(:,4)]); 
% axis([0 60 -1.5e-2  1.5e-2]);
xlabel('time [s]'); ylabel('rb/2V [deg]');
pbaspect([5 1 1]);

% figure(5)
subplot(5,1,5)
plot(t,[yt(:,5) yyt(:,5)]); 
% axis([0 60 -2 2]);
xlabel('time [s]'); ylabel('ay [m/s^2]');
pbaspect([5 1 1]);

pause
%% Spectral Analysis - Analytical

%simulation parameters
% w = logspace(-2,2,300);
fs = 1/dt;     % sample frequency
omega = 2*pi*fs*[0:2:(N/2)-1]/N;
omega2 = 2*pi*fs*[0:1:(N/2)-1]/N;

%controlled
temp = bode(Acon,B,Cexpanded(1,:),Dexpanded(1,:),4:5,omega);
temp = temp(:,1)+temp(:,2);
Sbeta  = temp.*temp;

temp = bode(Acon,B,Cexpanded(2,:),Dexpanded(2,:),4:5,omega); 
temp = temp(:,1)+temp(:,2);
Sphi   = temp.*temp;

temp = bode(Acon,B,Cexpanded(3,:),Dexpanded(3,:),4:5,omega); 
temp = temp(:,1)+temp(:,2);
Spp    = temp.*temp;

temp = bode(Acon,B,Cexpanded(4,:),Dexpanded(4,:),4:5,omega);
temp = temp(:,1)+temp(:,2);
Srr    = temp.*temp;

temp = bode(Acon,B,Cexpanded(5,:),Dexpanded(5,:),4:5,omega); 
temp = temp(:,1)+temp(:,2);
Say    = temp.*temp;

%reduced
temp = bode(Ared,Bred,Cred(1,:),Dred(1,:),4:5,omega);
temp = temp(:,1)+temp(:,2);
Sbetared  = temp.*temp;

temp = bode(Ared,Bred,Cred(2,:),Dred(2,:),4:5,omega); 
temp = temp(:,1)+temp(:,2);
Sphired   = temp.*temp;

temp = bode(Ared,Bred,Cred(3,:),Dred(3,:),4:5,omega); 
temp = temp(:,1)+temp(:,2);
Sppred    = temp.*temp;

temp = bode(Ared,Bred,Cred(4,:),Dred(4,:),4:5,omega);
temp = temp(:,1)+temp(:,2);
Srrred    = temp.*temp;

temp = bode(Ared,Bred,Cred(5,:),Dred(5,:),4:5,omega); 
temp = temp(:,1)+temp(:,2);
Sayred    = temp.*temp;

%synthesis
Sxx  = [Sbeta Sphi Spp Srr Say];
Sxxred = [Sbetared Sphired Sppred Srrred Sayred];


y     = lsim(Acon,B,Cexpanded,Dexpanded,w1,t) ... 
      + lsim(Acon,B,Cexpanded,Dexpanded,w2,t) ...
      + lsim(Acon,B,Cexpanded,Dexpanded,w3,t);
      
yred  = lsim(Ared,Bred,Cred,Dred,w1,t) ... 
      + lsim(Ared,Bred,Cred,Dred,w2,t) ...
      + lsim(Ared,Bred,Cred,Dred,w3,t);

%System responses
%controlled
beta  = y(:,1);
phi   = y(:,2);
p  = y(:,3);
r  = y(:,4);
ay    = y(:,5);

%reduced
betared  = yred(:,1);
phired   = yred(:,2);
pred     = yred(:,3);
rred     = yred(:,4);
ayred    = yred(:,5);

%% Spectral Analysis - Experimental - fft.m

%controlled
BETA  = fft(beta);
PHI   = fft(phi);
P     = fft(p);
R     = fft(r);
AY    = fft(ay);

Pbeta  = (dt/N)*( BETA.*conj(BETA));
Pphi   = (dt/N)*(  PHI.*conj(PHI));
Pp     = (dt/N)*(    P.*conj(P));
Pr     = (dt/N)*(    R.*conj(R));
Pay    = (dt/N)*(AY.*conj(AY));

Sfft = [Pbeta(1:round(N/2)-1)...
        Pphi(1:round(N/2)-1)...
        Pp(1:round(N/2)-1)...
        Pr(1:round(N/2)-1)...
        Pay(1:round(N/2)-1)];
     
     
%reduced
BETARED  = fft(betared);
PHIRED   = fft(phired);
PRED     = fft(pred);
RRED     = fft(rred);
AYRED    = fft(ayred);

Pbetared  = (dt/N)*( BETARED.*conj(BETARED));
Pphired   = (dt/N)*(  PHIRED.*conj(PHIRED));
Ppred     = (dt/N)*(    PRED.*conj(PRED));
Prred     = (dt/N)*(    RRED.*conj(RRED));
Payred    = (dt/N)*(AYRED.*conj(AYRED));

Sfftred = [Pbetared(1:round(N/2)-1)...
           Pphired(1:round(N/2)-1)...
           Ppred(1:round(N/2)-1)...
           Prred(1:round(N/2)-1)...
           Payred(1:round(N/2)-1)];

%% Spectral Analysis - Experimental - pwelch.m

%controlled
Sbetaw  = pwelch(beta,omega2,[],N,fs); Sbetaw = Sbetaw/2; Sbetaw = Sbetaw(1:end-1); 
Sphiw   = pwelch(phi,omega2,[],N,fs); Sphiw = Sphiw/2; Sphiw = Sphiw(1:end-1);
Spw     = pwelch(p,omega2,[],N,fs); Spw = Spw/2; Spw = Spw(1:end-1);
Srw     = pwelch(r,omega2,[],N,fs); Srw = Srw/2; Srw = Srw(1:end-1);
Sayw    = pwelch(ay,omega2,[],N,fs); Sayw = Sayw/2; Sayw = Sayw(1:end-1);

Sw = [Sbetaw Sphiw Spw Srw Sayw];


%reduced
Sbetaredw  = pwelch(betared,omega2,[],N,fs); Sbetaredw = Sbetaredw/2; Sbetaredw = Sbetaredw(1:end-1); 
Sphiredw   = pwelch(phired,omega2,[],N,fs); Sphiredw = Sphiredw/2; Sphiredw = Sphiredw(1:end-1);
Spredw     = pwelch(pred,omega2,[],N,fs); Spredw = Spredw/2; Spredw = Spredw(1:end-1);
Srredw     = pwelch(rred,omega2,[],N,fs); Srredw = Srredw/2; Srredw = Srredw(1:end-1);
Sayredw    = pwelch(ayred,omega2,[],N,fs); Sayredw = Sayredw/2; Sayredw = Sayredw(1:end-1);

Sredw = [Sbetaredw Sphiredw Spredw Srredw Sayredw];

%% Plotting

subplot(5,3,1)
loglog(omega,Sxx(:,1),'--',omega,Sxxred(:,1));
ylabel('Sbeta [deg^2*s]')
title('Analytical','FontSize',12)
legend('Controlled','Reduced')
grid on 
subplot(5,3,4)
loglog(omega,Sxx(:,2),'--',omega,Sxxred(:,2));
ylabel('Sphi [deg^2*s]')
grid on
subplot(5,3,7)
loglog(omega,Sxx(:,3),'--',omega,Sxxred(:,3));
ylabel('Spp [deg^2*s]')
grid on
subplot(5,3,10)
loglog(omega,Sxx(:,4),'--',omega,Sxxred(:,4));
ylabel('Srr [deg^2*s]')
grid on
subplot(5,3,13)
loglog(omega,Sxx(:,5),'--',omega,Sxxred(:,5));
xlabel('omega'); ylabel('Say [deg^2*s]')
grid on

subplot(5,3,2)
loglog(omega2,Sfft(:,1),'--',omega2,Sfftred(:,1));
title('Experimental - fft.m','FontSize',12)
grid on
subplot(5,3,5)
loglog(omega2,Sfft(:,2),'--',omega2,Sfftred(:,2));
grid on
subplot(5,3,8)
loglog(omega2,Sfft(:,3),'--',omega2,Sfftred(:,3));
grid on
subplot(5,3,11)
loglog(omega2,Sfft(:,4),'--',omega2,Sfftred(:,4));
grid on
subplot(5,3,14)
loglog(omega2,Sfft(:,5),'--',omega2,Sfftred(:,5));
xlabel('omega')
grid on

subplot(5,3,3)
loglog(omega2,Sw(:,1),'--',omega2,Sredw(:,1));
title('Experimental - pwelch.m','FontSize',12)
grid on
subplot(5,3,6)
loglog(omega2,Sw(:,2),'--',omega2,Sredw(:,2));
grid on
subplot(5,3,9)
loglog(omega2,Sw(:,3),'--',omega2,Sredw(:,3));
grid on
subplot(5,3,12)
loglog(omega2,Sw(:,4),'--',omega2,Sredw(:,4));
grid on
subplot(5,3,15)
loglog(omega2,Sw(:,5),'--',omega2,Sredw(:,5));
xlabel('omega')
grid on

pause

%% Variances

% Zero lists
Vara = zeros(1,5);
Varared = zeros(1,5);

Varfft = zeros(1,5);
Varfftred = zeros(1,5);

Varw = zeros(1,5);
Varwred = zeros(1,5);

%Integration
for i=1:length(omega)-1;
  for j=1:5
    Vara(j)=Vara(j)+(omega(i+1)-omega(i))*Sxx(i,j);
    Varared(j)=Varared(j)+(omega(i+1)-omega(i))*Sxxred(i,j);
    
    Varfft(j)=Varfft(j)+(omega(i+1)-omega(i))*Sfft(i,j);
    Varfftred(j)=Varfftred(j)+(omega(i+1)-omega(i))*Sfftred(i,j);
    
    Varw(j)=Varw(j)+(omega(i+1)-omega(i))*Sw(i,j);
    Varwred(j)=Varwred(j)+(omega(i+1)-omega(i))*Sredw(i,j);
  end
end

% Normalize 
Vara = Vara/pi
Varared = Varared/pi

Varfft = Varfft/pi
Varfftred = Varfftred/pi

Varw = Varw/pi
Varwred = Varwred/pi


% var.m method

% controlled
Vbeta = var(beta);
Vphi = var(phi);
Vp = var(p);
Vr = var(r);
Vay = var(ay);

Var = [Vbeta Vphi Vp Vr Vay]

% reduced
Vbetared = var(betared);
Vphired = var(phired);
Vpred = var(pred);
Vrred = var(rred);
Vayred = var(ayred);

Varred = [Vbetared Vphired Vpred Vrred Vayred]
