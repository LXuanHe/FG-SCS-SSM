%Main program
clc;clear;
K=0.9738;  %Scale factor
xlow=[2.5 3 6 50 11 0];
xupp=[5.5 12 30 130 39 3];
nd=size(xlow,2);

%% Sample points
A=dlmread('out50.txt');
A=unique(A,'rows','stable');
S=A(:,1:nd);           %Design variable
upTemI=A(1,nd+1);      %Inner side temperature rise
upTemO=A(1,nd+2);      %Outer side temperature rise
Y=A(:,nd+3);           %Collapse load

P=0;
psoiteration=0;     %Iterations
wholeiteration=0;
psoerror=1;         %Relative error
while psoiteration<20&&psoerror>=1e-2
    %% Optimization based on HSM and PSO
    hsmiteration=0; %Iterations
    hsmerror=1;     %Relative error
    while hsmiteration<=10&&hsmerror>=1e-2
        [xnew,ynew,Pco,mse,wt]=ego_kriging(S,Y,upTemI,upTemO,K,xlow,xupp);
        P=[P;Pco];
        %% Abaqus
        Preal=realfval(xnew,upTemI,upTemO);
        %% Relative error between HSM and abaqus
        hsmerror=abs(Preal-Pco)/Preal
        %% Update sample point
        S=[S;xnew];
        Y=[Y;Preal];
        hsmiteration=hsmiteration+1;
        wholeiteration=wholeiteration+1;
        %% Output
        EP=[Preal,hsmerror];
        C=[S,Y];
        dlmwrite('Error_Preal.txt',EP,'-append')
        dlmwrite('shuchu.txt',C)
        save('data')
    end
    %% PSO optimization
    dmodel=HSMmodel(S,Y,upTemI,upTemO,K,xlow,xupp);                   %Construct the better HSM
    [xopt,yopt,Fco,wt]=psofval(S,Y,upTemI,upTemO,K,xlow,xupp,dmodel); %Optimization based on PSO
    P=[P;Fco];
    Freal=realfval(xopt,upTemI,upTemO);                               %Calculation based on abaqus
    %% Relative error between HSM and abaqus
    psoerror=abs(Freal-Fco)/Freal
    psoiteration=psoiteration+1;
    %% Update sample point
    S=[S;xopt];
    Y=[Y;Freal];
    %% Output
    EF=[Freal,psoerror];
    C=[S,Y];
    dlmwrite('Error_Freal.txt',EF,'-append')
    dlmwrite('shuchu.txt',C)
    save('data')
end
save('ALLdata')