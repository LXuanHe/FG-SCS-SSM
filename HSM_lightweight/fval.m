function [c,ceq] = fval(x,upTemI,upTemO,K,dmodel,n)
%Optimization constraint

[Pimp,mse]=HSMfval(x,upTemI,upTemO,K,dmodel);

c=33014992.0-Pimp;
ceq=[];

end

