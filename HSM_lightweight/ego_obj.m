function [y]=ego_obj(x,upTemI,upTemO,K,dmodel,S)
%Objective function (Objective-pursuing learning method (OPLM))
[yr,mse]=HSMfval(x,upTemI,upTemO,K,dmodel);
yy=(abs(yr-33014992.0))/sqrt(mse);

wt=objval(x);
y=yy*wt^3;

end