%Main program
clc;clear;
%    (ts tr h Na Nc k)
xlow=[1 1.5 4 50 12 0];
xupp=[5 10 20 80 30 10];
n=size(xupp,2);  %Dimension
upTemI=400;      %Inner side temperature rise
upTemO=400;      %Outer side temperature rise
%% Run optimization based on the PSO algorithm
options = psooptimset('ConstrBoundary','soft','populationsize',40,'generations',150,'PlotFcns',@psoplotbestf);%,
[xopt,yopt]=pso(@(x) objval(x),n,[],[],[],[],xlow,xupp,@(x) fval(x,upTemI,upTemO),options);
xopt(n-2)=round(xopt(n-2));
xopt(n-1)=round(xopt(n-1));
%% SSM
Pcr=FGMSSM_T(xopt,upTemI,upTemO);
%% Abaqus
Preal=realfval(xopt,upTemI,upTemO);
%% Relative error
psoerror=abs(Preal-Pcr)/Preal
save('data')
