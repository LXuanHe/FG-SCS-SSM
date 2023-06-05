function [xnew,ynew,Pco,mse,wt]=ego_kriging(S,Y,upTemI,upTemO,K,xlow,xupp)
[m,n]=size(S);
[m,nc]=size(Y);
%% Construct HSM
dmodel=HSMmodel(S,Y,upTemI,upTemO,K,xlow,xupp);

 %% Run optimization based on the PSO algorithm
options = psooptimset('ConstrBoundary','soft','populationsize',100,'generations',200);%,'PlotFcns',@psoplotbestf
[xnew,ynew]=pso(@(x) ego_obj(x,upTemI,upTemO,K,dmodel,S),n,[],[],[],[],xlow,xupp,@(x) ego_constraint(x,upTemI,upTemO,K,dmodel),options);
xnew(n-2)=round(xnew(n-2));
xnew(n-1)=round(xnew(n-1));
%% The optimal solution is substituted into HSM to predict the ultimate load
[Pco,mse]=HSMfval(xnew,upTemI,upTemO,K,dmodel);
%% Calculate the optimal value(weight)
wt=objval(xnew);
end
