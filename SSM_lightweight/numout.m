function [F,M] = numout(wc)
%Data processing of abaqus results
if wc==1                              %Working condition
    f=textread('output1_F.txt','%s'); %Temperature load
    f=cell2mat(f(2,1));
    long_f=length(f);
    ff=f(:,1:long_f-1);
    F=str2num(ff);                    %Unit:N
elseif wc==2
    f=textread('output2_F.txt','%s'); %Critical buckling load
    f=cell2mat(f(5,1));
    long_f=length(f);
    ff=f(:,1:long_f);
    F=str2num(ff);                    %Unit:N
elseif wc==3
    f=textread('output3_F.txt','%s'); %Collapse load
    f=cell2mat(f);
    long_ff=length(f);
    ff=f(:,10:long_ff);
    F=str2num(ff);                    %Unit:N
end

rw=textread('output_W.txt','%s');     %Weight
ww=rw{1};
long_ww=length(ww);
jw=ww(:,4:long_ww);
M=str2num(jw);                        %Unit:kg
end

