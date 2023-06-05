function [c,ceq]=ego_constraint(x,upTemI,upTemO,K,dmodel)
%Constraint function
[yr,mse]=HSMfval(x,upTemI,upTemO,K,dmodel);  %HSM predicts collapse load

g=(33014992.0-yr);
c=g;
ceq=[];
end

