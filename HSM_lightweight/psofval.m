function [xopt,yopt,Fco,wt] = psofval(S,Y,upTemI,upTemO,K,xlow,xupp,dmodel)
%Optimization based on PSO 
[m,n]=size(S);
[m,nc]=size(Y);

%% Run optimization based on the PSO algorithm
options = psooptimset('ConstrBoundary','soft','populationsize',60,'generations',200);%,'PlotFcns',@psoplotbestf
[xopt,yopt]=pso(@(x) objval(x),n,[],[],[],[],xlow,xupp,@(x) fval(x,upTemI,upTemO,K,dmodel,n),options);
xopt(n-2)=round(xopt(n-2));
xopt(n-1)=round(xopt(n-1));

%% The optimal solution is substituted into HSM to predict the ultimate load
[Fco,mse]=HSMfval(xopt,upTemI,upTemO,K,dmodel);
%% Calculate the optimal value(weight)
wt=objval(xopt);

end