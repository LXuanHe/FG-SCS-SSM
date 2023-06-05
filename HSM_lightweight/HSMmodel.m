function [dmodel] = HSMmodel(S,Y,upTemI,upTemO,K,xlow,xupp)
%Construct the HSM
%The input and output of the sample are S and Y, respectively
[m,n]=size(S);
[m,nc]=size(Y);
for i=1:m
    g=FGMSSM_T(S(i,:),upTemI,upTemO);    %SSM
    error_g(i,1)=Y(i)-K*g;               %Errors between accurate models and SSM
end
dmodel=surrogate(S,error_g,1,xlow,xupp); %Construct the HSM
end

