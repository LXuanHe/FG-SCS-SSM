function [wt] = objval(x)
%Optimization objective
%The mass of the stiffened shell with variable section

% Design variables
ts=x(1);
tr=x(2);
h=x(3);
Na=round(x(4));
Nc=round(x(5));
k=x(6);

PI=acos(-1);
L=2000;  %Height
R=1500;  %Radius

miduc=5.7e-6;   %Inner side density,kg/mm^3
midum=4.429e-6; %Outer side density

syms z
V=(1/2+z/ts)^k;
midu=miduc*V+midum*(1-V);

w11=2*PI*(R+ts/2+z)*midu*L;
w1=integral(matlabFunction(w11),-ts/2,ts/2);   %Skin
w1=real(w1);
w2=Na*tr*h*L*miduc;                            %Axial stiffeners
w3=Nc*PI*(R^2-(R-h)^2)*tr*miduc;               %Circumferential stiffeners
wt=w1+w2+w3;

end

