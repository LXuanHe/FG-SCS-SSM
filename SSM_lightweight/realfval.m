function [ yreal ] = realfval(x,upTemI,upTemO)
%The critical buckling load calculated based on ABAQUS

x=[x,upTemI,upTemO];
dlmwrite('inputabaqus.txt',x')

system(['abaqus cae ','script',' noGUI','=E:\HSM_AK_FGMshell\SSM_lightweight\1.py']);
tf=numout(1);
dlmwrite('inputabaqus.txt',tf','-append') %Added temperature load
system(['abaqus cae ','script',' noGUI','=E:\HSM_AK_FGMshell\SSM_lightweight\2.py']);
lf=numout(2);

yreal=lf;

end

