function [c,ceq] = fval(x,upTemI,upTemO)
%Optimization constraint

Pcr=FGMSSM_T(x,upTemI,upTemO);

c=2496633-Pcr;
ceq=[];

end

