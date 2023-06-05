function [Pco,mse]=HSMfval(x,upTemI,upTemO,K,dmodel)
%HSM predicts collapse load

Pcr=FGMSSM_T(x,upTemI,upTemO); %SSM

ts=x(1);
tr=x(2);
h=x(3);
Na=round(x(4));
Nc=round(x(5));
k=x(6);

xx=[ts,tr,h,Na,Nc,k];
[g_preerror,dy,mse]=predictor(xx,dmodel); %Prediction error
Pco=K*Pcr+g_preerror;                     %Prediction collapse load

end