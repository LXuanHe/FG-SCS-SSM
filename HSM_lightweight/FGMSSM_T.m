function [Pcr]=FGMSSM_T(x,upTemI,upTemO)
%The critical buckling load of FG-SCS based on SSM

%% Design variables
ts=x(1);
tr=x(2);
h=x(3);
Na=round(x(4));
Nc=round(x(5));
k=x(6);

PI=acos(-1);
%% Structural physical parameters
L=2000.0;           %Height
R=1500.0;           %Radius

b=2*PI*R/Na;        %Axial spacing distance
d=L/(Nc-1);         %Circumferential spacing distance

mui=0.288;          %Inner side poisson ratio
muo=0.288;          %Outer side poisson ratio

%% Temperature
Tem=300;            %Initial
TemO=Tem+upTemO;    %Outer side
TemI=Tem+upTemI;    %Inner side
TemIO=TemI-TemO;

%% Temperature correlation coefficient
IE0=244.27e3;         %ZrO2,N/mm^2
IE=0;
IE1=-1.3707e-9;
IE2=1.21393e-12;
IE3=-3.6814e-16;

IA0=12.768e-6;        %1/K
IA=0;
IA1=-1.49e-3;
IA2=0.1e-5;
IA3=-0.6775e-11;

IK0=1.7;              %W/m*K
IK=0;
IK1=0.1276e-3;
IK2=0.66485e-5;
IK3=0;

OE0=122.56e3;         %Ti-6Al-4V,N/mm^2
OE=0;
OE1=-4.5866e-10;
OE2=0;
OE3=0;

OA0=7.5788e-6;        %1/K
OA=0;
OA1=6.5e-4;
OA2=0.31467e-6;
OA3=0;

OK0=1.20947;          %W/m*K
OK=0;
OK1=0.0139375;
OK2=0;
OK3=0;

%% Geometrical parameters of section
Ask=2*PI*R*ts;       %Skin cross section area
Ast_zhou=Na*tr*h;    %Axial stiffeners section area

Vsk=(2*b+tr)*(2*d+tr)*ts;                      %Skin volume of unit cell 
Vst=3*(2*b+tr)*tr*h+3*(2*d+tr)*tr*h-9*tr*tr*h; %Stiffeners volume of unit cell 

[Esk,Est,afsk,afst]=EPP(ts,k,TemI,TemO,TemIO); %Equivalent material parameter
E_sk=Esk;
E_st=3*tr*h*(b+d+tr)/Vst*Est;
h0=(ts+h)/2/(1+E_sk/E_st*Vsk/Vst)-ts/2;

%% Volume fraction
syms z
v=((z+h0)/ts+1);
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
Temshell=TemO+TemIO/C*(v-T1+T2-T3+T4-T5); %Skin temperature distribution

%% Material parameters at different temperatures
Ei=IE0*(IE/Temshell+1+IE1*Temshell+IE2*Temshell^2+IE3*Temshell^3);  %Inner side
afi=IA0*(IA/Temshell+1+IA1*Temshell+IA2*Temshell^2+IA3*Temshell^3);
Ki=IK0*(IK/Temshell+1+IK1*Temshell+IK2*Temshell^2+IK3*Temshell^3);

Eo=OE0*(OE/Temshell+1+OE1*Temshell+OE2*Temshell^2+OE3*Temshell^3);  %Outer side
afo=OA0*(OA/Temshell+1+OA1*Temshell+OA2*Temshell^2+OA3*Temshell^3);
Ko=OK0*(OK/Temshell+1+OK1*Temshell+OK2*Temshell^2+OK3*Temshell^3);

Eskin=Ei*V+Eo*(1-V);
afskin=afi*V+afo*(1-V);
Kskin=Ki*V+Ko*(1-V);

%% Temperature load
if TemO~=TemI
    dLsk=afskin*(Temshell-Tem);
    dLsk=integral(matlabFunction(dLsk),-(ts+h0),-h0)/ts;
else
    dLsk=afskin*upTemO;
    dLsk=integral(matlabFunction(dLsk),-(ts+h0),-h0)/ts;
end
dLst=afst*upTemI;
dL=dLsk-(dLsk-dLst)/(1+Ask*Esk/Ast_zhou/Est);
Ftem=Esk*dL*Ask+Est*dL*Ast_zhou;

%% Skin stiffness
A11skin=Eskin/(1-mui*mui);
A11skin=integral(matlabFunction(A11skin),-(ts+h0),-h0);
A22skin=A11skin;
A12skin=Eskin*mui/(1-mui*mui);
A12skin=integral(matlabFunction(A12skin),-(ts+h0),-h0);
A21skin=A12skin;
A66skin=Eskin/2/(1+mui);
A66skin=integral(matlabFunction(A66skin),-(ts+h0),-h0);

B11skin=Eskin*z/(1-mui*mui);
B11skin=integral(matlabFunction(B11skin),-(ts+h0),-h0);
B22skin=B11skin;
B12skin=Eskin*z*mui/(1-mui*mui);
B12skin=integral(matlabFunction(B12skin),-(ts+h0),-h0);
B21skin=B12skin;
B66skin=Eskin*z/2/(1+mui);
B66skin=integral(matlabFunction(B66skin),-(ts+h0),-h0);

D11skin=Eskin*z*z/(1-mui*mui);
D11skin=integral(matlabFunction(D11skin),-(ts+h0),-h0);
D22skin=D11skin;
D12skin=Eskin*z*z*mui/(1-mui*mui);
D12skin=integral(matlabFunction(D12skin),-(ts+h0),-h0);
D21skin=D12skin;
D66skin=Eskin*z*z/2/(1+mui);
D66skin=integral(matlabFunction(D66skin),-(ts+h0),-h0);

%% Stiffeners stiffness
Gs=Est/2.0/(1+mui);
Gr=Gs;
As=tr*h;
Ar=As;
Zs=h/2-h0;
Zr=Zs;

Is=tr*h^3/12;
Ir=Is;
Js=h*tr^3/16*(16/3-3.36*tr/h*(1-tr^4/12/h^4));
Jr=Js;

A11rib=Est*As/b;
A22rib=Est*Ar/d;
A12rib=0;
A21rib=0;
A66rib=0;
a66rib=1/12*b*b*d/Est/Ir+1/12*b*d*d/Est/Is+d/Gs/As/(6/5)+b/Gs/Ar/(6/5);
A66rib=1/a66rib;

B11rib=Zs*Est*As/b;
B22rib=Zr*Est*Ar/d;
B12rib=0;
B21rib=0;
B66rib=0;

D11rib=Est*(Is+Zs*Zs*As)/b;
D22rib=Est*(Ir+Zr*Zr*Ar)/d;
D12rib=0;
D21rib=0;
D66rib=Gs*Js*(1/b+1/d)/2;

%% Local buckling load
D=Esk*ts^3/12/(1-mui^2);
mm1=10000000;
nn1=10000000; 
N_plate=1e15;
kp=(ts*b+tr*h)/ts/b;
for m=1:99
    for n=1:99
        Nx=PI^2*d^2*D*(m^2/d^2+n^2/b^2)^2/m^2*kp;
        if Nx<N_plate
            N_plate=Nx;
            mm1=m;
            nn1=n;
        end
    end
end
Pcr_local=2*PI*R*N_plate-Ftem;

%% Global critical buckling load
A11=A11skin+A11rib;
A22=A22skin+A22rib;
A12=A12skin+A12rib;
A21=A21skin+A21rib;
A66=A66skin+A66rib;

B11=B11skin+B11rib;
B22=B22skin+B22rib;
B12=B12skin+B12rib;
B21=B21skin+B21rib;
B66=B66skin+B66rib;

D11=D11skin+D11rib;
D22=D22skin+D22rib;
D12=D12skin+D12rib;
D21=D21skin+D21rib;
D66=D66skin+D66rib;

ma=1e10;
nb=1e10;
mm2=10000000;
nn2=10000000;
N_all=1e15;
for m=1:99
    for n=1:99
        ma=m*PI/L;
        nb=n/R;
        T11=A11*ma^2+A66*nb^2;
        T22=A22*nb^2+A66*ma^2;
        T33=A22/R^2-(2*B12*ma^2+2*B22*nb^2)/R+D11*ma^4+D22*nb^4+(2*D12+4*D66)*ma^2*nb^2;
        T12=(A12+A66)*ma*nb;
        T21=T12;
        T13=A12/R*ma-B11*ma^3-(B12+2*B66)*ma*nb^2;
        T31=T13;
        T23=A22/R*nb-B22*nb^3-(B12+2*B66)*ma^2*nb;
        T32=T23;
        
        TU=[T11 T12 T13;T21 T22 T23;T31 T32 T33];
        TD=[T11 T12;T21 T22];
        
        Nx=det(TU)/det(TD)/ma/ma;
        if Nx<N_all
            N_all=Nx;
            mm2=m;
            nn2=n;
        end
    end
end

Pcr_all=2*PI*R*N_all-Ftem;

% fprintf('Local buckling load is %g\n %g\n %g\n',Pcr_local,mm1,nn1);
% fprintf('Global critical buckling load is %g\n %g\n %g\n',Pcr_all,mm2,nn2);
% fprintf('Temperature load is %g\n',Ftem);

Pcr=min(Pcr_all,Pcr_local);
% Pcr=Pcr_all;
end
    
