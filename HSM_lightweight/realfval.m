function [ yreal ] = realfval(x,upTemI,upTemO)
%The critical buckling load calculated based on ABAQUS
x=[x,upTemI,upTemO];
dlmwrite('inputabaqus.txt',x')
tic
system(['abaqus cae ','script',' noGUI','=E:\HSM_AK_FGMshell\HSM_lightweight\1.py']);
tf=numout(1);
dlmwrite('inputabaqus.txt',tf','-append') %Added temperature load
system(['abaqus cae ','script',' noGUI','=E:\HSM_AK_FGMshell\HSM_lightweight\3.py']);
jf=numout(3);
toc
yreal=jf;

end

