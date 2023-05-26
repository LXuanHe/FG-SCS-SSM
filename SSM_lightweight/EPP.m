function [Esk,Est,afsk,afst] = EPP(ts,k,TemI,TemO,TemIO)
%Equivalent physical proerty

OE0=349.55e3;         %Al2O3,N/mm^2
OE=0;
OE1=-3.853e-10;
OE2=4.027e-13;
OE3=-1.673e-16;

OA0=6.8269e-6;        %1/K
OA=0;
OA1=1.838e-4;
OA2=0;
OA3=0;

OK0=-14.087;          %W/m*K
OK=-1123.6;
OK1=-6.227e-3;
OK2=0;
OK3=0;

IE0=122.56e3;         %Ti-6Al-4V,N/mm^2
IE=0;
IE1=-4.5866e-10;
IE2=0;
IE3=0;

IA0=7.5788e-6;        %1/K
IA=0;
IA1=6.5e-4;
IA2=0.31467e-6;
IA3=0;

IK0=1.20947;          %W/m*K
IK=0;
IK1=0.0139375;
IK2=0;
IK3=0;

%% Volume fraction
syms z
v=(z/ts+1/2);
V=v^k;

%% Distribution of temperature
KI=IK0*(IK/TemI+1+IK1*TemI+IK2*TemI^2+IK3*TemI^3);
KO=OK0*(OK/TemO+1+OK1*TemO+OK2*TemO^2+OK3*TemO^3);
KIO=KI-KO;
C1=KIO/(k+1)/KO;
C2=KIO^2/(2*k+1)/KO^2;
C3=KIO^3/(3*k+1)/KO^3;
C4=KIO^4/(4*k+1)/KO^4;
C5=KIO^5/(5*k+1)/KO^5;
C=1-C1+C2-C3+C4-C5;
T1=KIO/(k+1)/KO*v^(k+1);
T2=KIO^2/(2*k+1)/KO^2*v^(2*k+1);
T3=KIO^3/(3*k+1)/KO^3*v^(3*k+1);
T4=KIO^4/(4*k+1)/KO^4*v^(4*k+1);
T5=KIO^5/(5*k+1)/KO^5*v^(5*k+1);
Temshell=TemO+TemIO/C*(v-T1+T2-T3+T4-T5);

%% Material parameters at different temperatures
Ei=IE0*(IE/Temshell+1+IE1*Temshell+IE2*Temshell^2+IE3*Temshell^3);
afi=IA0*(IA/Temshell+1+IA1*Temshell+IA2*Temshell^2+IA3*Temshell^3);
Ki=IK0*(IK/Temshell+1+IK1*Temshell+IK2*Temshell^2+IK3*Temshell^3);

Eo=OE0*(OE/Temshell+1+OE1*Temshell+OE2*Temshell^2+OE3*Temshell^3);
afo=OA0*(OA/Temshell+1+OA1*Temshell+OA2*Temshell^2+OA3*Temshell^3);
Ko=OK0*(OK/Temshell+1+OK1*Temshell+OK2*Temshell^2+OK3*Temshell^3);

Eskin=Ei*V+Eo*(1-V);
afskin=afi*V+afo*(1-V);
Kskin=Ki*V+Ko*(1-V);

%Skin
Esk=integral(matlabFunction(Eskin),-ts/2,ts/2)/ts;
afsk=integral(matlabFunction(afskin),-ts/2,ts/2)/ts;
%Stiffeners
Est=IE0*(IE/TemI+1+IE1*TemI+IE2*TemI^2+IE3*TemI^3);
afst=IA0*(IA/TemI+1+IA1*TemI+IA2*TemI^2+IA3*TemI^3);

end

