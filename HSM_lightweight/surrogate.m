function dmodel=surrogate(S,Y,type,xmin,xmax)
%Construct surrogate model
%S and Y are inputs and outputs
%Type equal to 1 is kriging
%xmin and xmax are upper and lower bounds on the input
[m,n]=size(S);
if type==1
    theta=10*ones(1,n);
    lob=0.001*ones(1,n);
    upb=20*ones(1,n);
    [dmodel, perf] = ...
        dacefit(S, Y, @regpoly1, @corrgauss, theta, lob, upb);
end

end
    

