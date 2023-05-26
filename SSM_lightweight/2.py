#!/usr/bin/env python
#coding=utf-8

#FGM����ֵ����

from abaqus import*
from part import*
from material import*
from section import*
from assembly import*
from step import*
from interaction import*
from load import*
from mesh import*
from job import*
from sketch import*
from visualization import*
from string import*
from math import*
from collections import OrderedDict

import section 
import regionToolset 
import displayGroupOdbToolset as dgo
import __main__
global PI
PI=acos(-1)

######## 
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')   #r��ʾ��ֻ����ʽ���ļ�  
HOU=float(fd.readlines()[0])
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')
RibWidth=float(fd.readlines()[1])
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')   
RibHeight=float(fd.readlines()[2])                                        #�����߶�
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')
VerticalRibNum=int(fd.readlines()[3])                                     #�ݽ���  90
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')
HorizontalRibNum=int(fd.readlines()[4])                                   #������  25
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')
vk=float(fd.readlines()[5])                                               #�������ָ��
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')
upTemI=float(fd.readlines()[6])                                           #�ڱ���������K
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')
upTemO=float(fd.readlines()[7])                                           #�����������K
fd.close()
fd=file('E:/HSM_AK_FGMshell/SSM_lightweight/inputabaqus.txt','r')
Jzhz=float(fd.readlines()[8])                                             #��Ч�¶Ⱥ��أ���λ:N
fd.close()

msize1=10               #����ߴ�
msize2=2.5              #�ߴ�2�����˼���ʱ��
computTime=0.2

Tem=300.0               #ʼ��
TemO=Tem+upTemO         #����棨�մɲ㣩����
TemI=Tem+upTemI         #�ڱ��棨�����㣩����
TemIO=TemI-TemO         #������������²�

cengshu=20              #�ֲ���� 
CH=(1.0/cengshu)*HOU    #������
InPt=3                  #������ֵ���
TemPt=1+cengshu*(InPt-1)#���¶ȵ���

GAO=750.0    
BANJING=500.0

global R
R=float(BANJING)
global H
H=float(GAO)
global gg


vn=VerticalRibNum        #�ݽ���
hn=HorizontalRibNum      #������
tw=RibWidth              #�������
h=RibHeight              #�����߶�
L=2*PI*R                 #�����ܳ�
vgap=L/vn                #�ݽ��
hgap=H/(hn-1)            #�����

RibAngle=30              #δ�õ�
alpha=RibAngle*pi/180.0  #30�㣬δ�õ�
ribtype=HorizontalRibNum #���ƻ���������Ի��������������ǽṹ��δ�õ�

ptname='p1'
modelname='aa'
computmethod='Linear_Buck'


##########�������#############
posongb=0.288
Iposongb=posongb     #�ڱ��棨���������ɱ�
Oposongb=posongb     #����棨�մɣ����ɱ�

Imidu=4.429e-9       #�ڱ��棨�������ܶȣ�g/mm^3
Omidu=3.97e-9        #����棨�մɣ��ܶ�

###############�¶�ϵ��###################
endTem_300=300.0     #Ԥ����������¶ȣ�K
endTem_325=325.0
endTem_350=350.0
endTem_375=375.0
endTem_400=400.0
endTem_425=425.0
endTem_450=450.0
endTem_475=475.0
endTem_500=500.0
endTem_525=525.0
endTem_550=550.0
endTem_575=575.0
endTem_600=600.0
endTem_625=625.0
endTem_650=650.0
endTem_675=675.0
endTem_700=700.0
endTem_725=725.0
endTem_750=750.0
endTem_775=775.0
endTem_800=800.0
endTem_825=825.0
endTem_850=850.0
endTem_875=875.0
endTem_900=900.0

OE0=349.55e3         #����棨�մɣ���ģ,��λ��N/mm^2 = MPa
OE=0.0
OE1=-3.853e-10
OE2=4.027e-13
OE3=-1.673e-16

OA0=6.8269e-6        #����棨�մɣ����ͣ���λ��1/K
OA=0.0
OA1=1.838e-4
OA2=0.0
OA3=0.0

OK0=-14.087          #����棨�մɣ��ȵ��ʣ���λ��W/m*K����/�ס���
OK=-1123.6
OK1=-6.227e-3
OK2=0.0
OK3=0.0

IE0=122.56e3         #�ڱ��棨��������ģ
IE=0.0
IE1=-4.5866e-10
IE2=0.0
IE3=0.0

IA0=7.5788e-6        #�ڱ��棨����������
IA=0.0
IA1=6.5e-4
IA2=0.31467e-6
IA3=0.0

IK0=1.20947          #�ڱ��棨�������ȵ���
IK=0.0
IK1=0.0139375
IK2=0.0
IK3=0.0

MY0=1540.9254        #������������,MPa
MY=0.0
MY1=-7.2038e-4
MY2=1.0429e-7
MY3=0.0

MH0=2306.1436        #����ǿ��ģ��,MPa
MH=0.0
MH1=-4.9993e-4
MH2=0.0
MH3=0.0


##########��Ƥ�¶�#############
Ki=IK0*(IK/TemI+1.0+IK1*TemI+IK2*TemI**2+IK3*TemI**3)
Ko=OK0*(OK/TemO+1.0+OK1*TemO+OK2*TemO**2+OK3*TemO**3)
Kio=Ki-Ko
C1=Kio/(vk+1)/Ko
C2=Kio**2/(2*vk+1)/Ko**2
C3=Kio**3/(3*vk+1)/Ko**3
C4=Kio**4/(4*vk+1)/Ko**4
C5=Kio**5/(5*vk+1)/Ko**5
C=-C1+C2-C3+C4-C5+1.0


zz=[]
for i in range(1, TemPt+1):                                 #ע�⣬for i in range(1,3) ,������1��2����û3
    zz.append(-HOU/2.0+(int(i)-1)/(float(TemPt-1))*HOU)     #������浽�ڱ���


Vv=[]
for i in range(0, TemPt): 
    Vv.append((zz[i]/HOU+0.5))


T1=[]
for i in range(0, TemPt): 
    T1.append(Kio/(vk+1.0)/Ko*Vv[i]**(vk+1.0))


T2=[]
for i in range(0, TemPt): 
    T2.append(Kio**2.0/(2.0*vk+1.0)/Ko**2*Vv[i]**(2.0*vk+1.0))


T3=[]
for i in range(0, TemPt): 
    T3.append(Kio**3.0/(3.0*vk+1.0)/Ko**3.0*Vv[i]**(3.0*vk+1.0))


T4=[]
for i in range(0, TemPt): 
    T4.append(Kio**4.0/(4.0*vk+1.0)/Ko**4.0*Vv[i]**(4.0*vk+1.0))


T5=[]
for i in range(0, TemPt): 
    T5.append(Kio**5.0/(5.0*vk+1.0)/Ko**5.0*Vv[i]**(5.0*vk+1.0))


Temshell=[]
for i in range(0, TemPt): 
    Temshell.append(TemO+TemIO/C*(Vv[i]-T1[i]+T2[i]-T3[i]+T4[i]-T5[i]))   #��Ƥ����TemPt���ο��¶ȵ㣬������浽�ڱ���Ⱦ�ֲ�

Temlist=[]
for i in range(1, TemPt+1): 
    Temlist.append(Temshell[TemPt-i])  #���ڱ��浽�����Ⱦ�ֲ��������ں���ʩ��ģ�飬ָ����Ƥ�¶ȷֲ��Ǵ��µ��浽�ϱ��棬���ڱڵ���ڣ�
                                       #�����Ҫ���¶ȵ������Ų����з��䣬���¶ȵ���������������ģ�����Ӱ���С����

#############�¶�300K�µĲ����������#############
Ei_300=IE0*(IE/endTem_300+1.0+IE1*endTem_300+IE2*endTem_300**2+IE3*endTem_300**3)
Eo_300=OE0*(OE/endTem_300+1.0+OE1*endTem_300+OE2*endTem_300**2+OE3*endTem_300**3)
alphai_300=IA0*(IA/endTem_300+1.0+IA1*endTem_300+IA2*endTem_300**2+IA3*endTem_300**3)
alphao_300=OA0*(OA/endTem_300+1.0+OA1*endTem_300+OA2*endTem_300**2+OA3*endTem_300**3)
Ki_300=IK0*(IK/endTem_300+1.0+IK1*endTem_300+IK2*endTem_300**2+IK3*endTem_300**3)
Ko_300=OK0*(OK/endTem_300+1.0+OK1*endTem_300+OK2*endTem_300**2+OK3*endTem_300**3)
YS_300=MY0*(MY/endTem_300+1.0+MY1*endTem_300+MY2*endTem_300**2+MY3*endTem_300**3)
HS_300=MH0*(MH/endTem_300+1.0+MH1*endTem_300+MH2*endTem_300**2+MH3*endTem_300**3)

zzz=[]
for i in range(1, cengshu+1):
    zzz.append(HOU/2.0+(-int(i)+0.5)/cengshu*HOU)   #ע�⣬for i in range(1,3) ,������1��2����û3
                                                    #��z���򣨺�ȣ�λ�ú��������ڵ���

Volumefenshu=[]   #�ڱ�����ϣ��մɣ���������������ڲ㵽��㣨��1��0��
for i in range(0, cengshu): 
    Volumefenshu.append((zzz[i]/HOU+0.5)**vk)


midu=[]           #�ܶȣ����¶��޹ء�������йأ�g/mm^3
for i in range(0, cengshu):         
    midu.append(Imidu*Volumefenshu[i]+Omidu*(1.0+(-Volumefenshu[i])))    #���ڣ��մɣ����ɵ��⣨������,���浯ģ�����ͣ���������һ��


tanxine_300=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_300.append(Ei_300*Volumefenshu[i]+Eo_300*(1.0+(-Volumefenshu[i])))


pengzhang_300=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_300.append(alphai_300*Volumefenshu[i]+alphao_300*(1.0+(-Volumefenshu[i])))
    

daore_300=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_300.append(Ki_300*Volumefenshu[i]+Ko_300*(1.0+(-Volumefenshu[i])))

   
qufu_300_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_300_a.append(YS_300*(1.0+(-Volumefenshu[i])+Ei_300/Eo_300*Volumefenshu[i]))

    
qianghua_300=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_300.append(Ei_300*Volumefenshu[i]+HS_300*(1.0+(-Volumefenshu[i])))


qufu_300_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_300_b.append(qufu_300_a[i]+0.1*qianghua_300[i])
    
    
YS_300_a=YS_300   #Ti-6Al-4V
YS_300_b=YS_300_a+0.1*HS_300


#############�¶�325K�µĲ����������###############
Ei_325=IE0*(IE/endTem_325+1.0+IE1*endTem_325+IE2*endTem_325**2+IE3*endTem_325**3)
Eo_325=OE0*(OE/endTem_325+1.0+OE1*endTem_325+OE2*endTem_325**2+OE3*endTem_325**3)
alphai_325=IA0*(IA/endTem_325+1.0+IA1*endTem_325+IA2*endTem_325**2+IA3*endTem_325**3)
alphao_325=OA0*(OA/endTem_325+1.0+OA1*endTem_325+OA2*endTem_325**2+OA3*endTem_325**3)
Ki_325=IK0*(IK/endTem_325+1.0+IK1*endTem_325+IK2*endTem_325**2+IK3*endTem_325**3)
Ko_325=OK0*(OK/endTem_325+1.0+OK1*endTem_325+OK2*endTem_325**2+OK3*endTem_325**3)
YS_325=MY0*(MY/endTem_325+1.0+MY1*endTem_325+MY2*endTem_325**2+MY3*endTem_325**3)
HS_325=MH0*(MH/endTem_325+1.0+MH1*endTem_325+MH2*endTem_325**2+MH3*endTem_325**3)

tanxine_325=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_325.append(Ei_325*Volumefenshu[i]+Eo_325*(1.0+(-Volumefenshu[i])))


pengzhang_325=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_325.append(alphai_325*Volumefenshu[i]+alphao_325*(1.0+(-Volumefenshu[i])))
    
    
daore_325=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_325.append(Ki_325*Volumefenshu[i]+Ko_325*(1.0+(-Volumefenshu[i])))

    
qufu_325_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_325_a.append(YS_325*(1.0+(-Volumefenshu[i])+Ei_325/Eo_325*Volumefenshu[i]))

    
qianghua_325=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_325.append(Ei_325*Volumefenshu[i]+HS_325*(1.0+(-Volumefenshu[i])))

    
qufu_325_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_325_b.append(qufu_325_a[i]+0.1*qianghua_325[i])
    

YS_325_a=YS_325   #Ti-6Al-4V
YS_325_b=YS_325_a+0.1*HS_325


#############�¶�350K�µĲ����������###############
Ei_350=IE0*(IE/endTem_350+1.0+IE1*endTem_350+IE2*endTem_350**2+IE3*endTem_350**3)
Eo_350=OE0*(OE/endTem_350+1.0+OE1*endTem_350+OE2*endTem_350**2+OE3*endTem_350**3)
alphai_350=IA0*(IA/endTem_350+1.0+IA1*endTem_350+IA2*endTem_350**2+IA3*endTem_350**3)
alphao_350=OA0*(OA/endTem_350+1.0+OA1*endTem_350+OA2*endTem_350**2+OA3*endTem_350**3)
Ki_350=IK0*(IK/endTem_350+1.0+IK1*endTem_350+IK2*endTem_350**2+IK3*endTem_350**3)
Ko_350=OK0*(OK/endTem_350+1.0+OK1*endTem_350+OK2*endTem_350**2+OK3*endTem_350**3)
YS_350=MY0*(MY/endTem_350+1.0+MY1*endTem_350+MY2*endTem_350**2+MY3*endTem_350**3)
HS_350=MH0*(MH/endTem_350+1.0+MH1*endTem_350+MH2*endTem_350**2+MH3*endTem_350**3)

tanxine_350=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_350.append(Ei_350*Volumefenshu[i]+Eo_350*(1.0+(-Volumefenshu[i])))


pengzhang_350=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_350.append(alphai_350*Volumefenshu[i]+alphao_350*(1.0+(-Volumefenshu[i])))
    
    
daore_350=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_350.append(Ki_350*Volumefenshu[i]+Ko_350*(1.0+(-Volumefenshu[i])))

    
qufu_350_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_350_a.append(YS_350*(1.0+(-Volumefenshu[i])+Ei_350/Eo_350*Volumefenshu[i]))

    
qianghua_350=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_350.append(Ei_350*Volumefenshu[i]+HS_350*(1.0+(-Volumefenshu[i])))

   
qufu_350_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_350_b.append(qufu_350_a[i]+0.1*qianghua_350[i])
    

YS_350_a=YS_350   #Ti-6Al-4V
YS_350_b=YS_350_a+0.1*HS_350

   
#############�¶�375K�µĲ����������###############
Ei_375=IE0*(IE/endTem_375+1.0+IE1*endTem_375+IE2*endTem_375**2+IE3*endTem_375**3)
Eo_375=OE0*(OE/endTem_375+1.0+OE1*endTem_375+OE2*endTem_375**2+OE3*endTem_375**3)
alphai_375=IA0*(IA/endTem_375+1.0+IA1*endTem_375+IA2*endTem_375**2+IA3*endTem_375**3)
alphao_375=OA0*(OA/endTem_375+1.0+OA1*endTem_375+OA2*endTem_375**2+OA3*endTem_375**3)
Ki_375=IK0*(IK/endTem_375+1.0+IK1*endTem_375+IK2*endTem_375**2+IK3*endTem_375**3)
Ko_375=OK0*(OK/endTem_375+1.0+OK1*endTem_375+OK2*endTem_375**2+OK3*endTem_375**3)
YS_375=MY0*(MY/endTem_375+1.0+MY1*endTem_375+MY2*endTem_375**2+MY3*endTem_375**3)
HS_375=MH0*(MH/endTem_375+1.0+MH1*endTem_375+MH2*endTem_375**2+MH3*endTem_375**3)

tanxine_375=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_375.append(Ei_375*Volumefenshu[i]+Eo_375*(1.0+(-Volumefenshu[i])))


pengzhang_375=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_375.append(alphai_375*Volumefenshu[i]+alphao_375*(1.0+(-Volumefenshu[i])))
    
    
daore_375=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_375.append(Ki_375*Volumefenshu[i]+Ko_375*(1.0+(-Volumefenshu[i])))

   
qufu_375_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_375_a.append(YS_375*(1.0+(-Volumefenshu[i])+Ei_375/Eo_375*Volumefenshu[i]))

    
qianghua_375=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_375.append(Ei_375*Volumefenshu[i]+HS_375*(1.0+(-Volumefenshu[i])))
 
  
qufu_375_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_375_b.append(qufu_375_a[i]+0.1*qianghua_375[i])


YS_375_a=YS_375   #Ti-6Al-4V
YS_375_b=YS_375_a+0.1*HS_375


#############�¶�400K�µĲ����������###############
Ei_400=IE0*(IE/endTem_400+1.0+IE1*endTem_400+IE2*endTem_400**2+IE3*endTem_400**3)
Eo_400=OE0*(OE/endTem_400+1.0+OE1*endTem_400+OE2*endTem_400**2+OE3*endTem_400**3)
alphai_400=IA0*(IA/endTem_400+1.0+IA1*endTem_400+IA2*endTem_400**2+IA3*endTem_400**3)
alphao_400=OA0*(OA/endTem_400+1.0+OA1*endTem_400+OA2*endTem_400**2+OA3*endTem_400**3)
Ki_400=IK0*(IK/endTem_400+1.0+IK1*endTem_400+IK2*endTem_400**2+IK3*endTem_400**3)
Ko_400=OK0*(OK/endTem_400+1.0+OK1*endTem_400+OK2*endTem_400**2+OK3*endTem_400**3)
YS_400=MY0*(MY/endTem_400+1.0+MY1*endTem_400+MY2*endTem_400**2+MY3*endTem_400**3)
HS_400=MH0*(MH/endTem_400+1.0+MH1*endTem_400+MH2*endTem_400**2+MH3*endTem_400**3)

tanxine_400=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_400.append(Ei_400*Volumefenshu[i]+Eo_400*(1.0+(-Volumefenshu[i])))


pengzhang_400=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_400.append(alphai_400*Volumefenshu[i]+alphao_400*(1.0+(-Volumefenshu[i])))
    
    
daore_400=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_400.append(Ki_400*Volumefenshu[i]+Ko_400*(1.0+(-Volumefenshu[i])))

    
qufu_400_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_400_a.append(YS_400*(1.0+(-Volumefenshu[i])+Ei_400/Eo_400*Volumefenshu[i]))

    
qianghua_400=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_400.append(Ei_400*Volumefenshu[i]+HS_400*(1.0+(-Volumefenshu[i])))

   
qufu_400_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_400_b.append(qufu_400_a[i]+0.1*qianghua_400[i])


YS_400_a=YS_400   #Ti-6Al-4V
YS_400_b=YS_400_a+0.1*HS_400


#############�¶�425K�µĲ����������###############
Ei_425=IE0*(IE/endTem_425+1.0+IE1*endTem_425+IE2*endTem_425**2+IE3*endTem_425**3)
Eo_425=OE0*(OE/endTem_425+1.0+OE1*endTem_425+OE2*endTem_425**2+OE3*endTem_425**3)
alphai_425=IA0*(IA/endTem_425+1.0+IA1*endTem_425+IA2*endTem_425**2+IA3*endTem_425**3)
alphao_425=OA0*(OA/endTem_425+1.0+OA1*endTem_425+OA2*endTem_425**2+OA3*endTem_425**3)
Ki_425=IK0*(IK/endTem_425+1.0+IK1*endTem_425+IK2*endTem_425**2+IK3*endTem_425**3)
Ko_425=OK0*(OK/endTem_425+1.0+OK1*endTem_425+OK2*endTem_425**2+OK3*endTem_425**3)
YS_425=MY0*(MY/endTem_425+1.0+MY1*endTem_425+MY2*endTem_425**2+MY3*endTem_425**3)
HS_425=MH0*(MH/endTem_425+1.0+MH1*endTem_425+MH2*endTem_425**2+MH3*endTem_425**3)

tanxine_425=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_425.append(Ei_425*Volumefenshu[i]+Eo_425*(1.0+(-Volumefenshu[i])))


pengzhang_425=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_425.append(alphai_425*Volumefenshu[i]+alphao_425*(1.0+(-Volumefenshu[i])))
    
    
daore_425=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_425.append(Ki_425*Volumefenshu[i]+Ko_425*(1.0+(-Volumefenshu[i])))

   
qufu_425_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_425_a.append(YS_425*(1.0+(-Volumefenshu[i])+Ei_425/Eo_425*Volumefenshu[i]))

   
qianghua_425=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_425.append(Ei_425*Volumefenshu[i]+HS_425*(1.0+(-Volumefenshu[i])))


qufu_425_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_425_b.append(qufu_425_a[i]+0.1*qianghua_425[i])
    

YS_425_a=YS_425   #Ti-6Al-4V
YS_425_b=YS_425_a+0.1*HS_425


#############�¶�450K�µĲ����������###############
Ei_450=IE0*(IE/endTem_450+1.0+IE1*endTem_450+IE2*endTem_450**2+IE3*endTem_450**3)
Eo_450=OE0*(OE/endTem_450+1.0+OE1*endTem_450+OE2*endTem_450**2+OE3*endTem_450**3)
alphai_450=IA0*(IA/endTem_450+1.0+IA1*endTem_450+IA2*endTem_450**2+IA3*endTem_450**3)
alphao_450=OA0*(OA/endTem_450+1.0+OA1*endTem_450+OA2*endTem_450**2+OA3*endTem_450**3)
Ki_450=IK0*(IK/endTem_450+1.0+IK1*endTem_450+IK2*endTem_450**2+IK3*endTem_450**3)
Ko_450=OK0*(OK/endTem_450+1.0+OK1*endTem_450+OK2*endTem_450**2+OK3*endTem_450**3)
YS_450=MY0*(MY/endTem_450+1.0+MY1*endTem_450+MY2*endTem_450**2+MY3*endTem_450**3)
HS_450=MH0*(MH/endTem_450+1.0+MH1*endTem_450+MH2*endTem_450**2+MH3*endTem_450**3)

tanxine_450=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_450.append(Ei_450*Volumefenshu[i]+Eo_450*(1.0+(-Volumefenshu[i])))


pengzhang_450=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_450.append(alphai_450*Volumefenshu[i]+alphao_450*(1.0+(-Volumefenshu[i])))
    
    
daore_450=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_450.append(Ki_450*Volumefenshu[i]+Ko_450*(1.0+(-Volumefenshu[i])))

   
qufu_450_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_450_a.append(YS_450*(1.0+(-Volumefenshu[i])+Ei_450/Eo_450*Volumefenshu[i]))

   
qianghua_450=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_450.append(Ei_450*Volumefenshu[i]+HS_450*(1.0+(-Volumefenshu[i])))


qufu_450_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_450_b.append(qufu_450_a[i]+0.1*qianghua_450[i])


YS_450_a=YS_450   #Ti-6Al-4V
YS_450_b=YS_450_a+0.1*HS_450


#############�¶�475K�µĲ����������###############
Ei_475=IE0*(IE/endTem_475+1.0+IE1*endTem_475+IE2*endTem_475**2+IE3*endTem_475**3)
Eo_475=OE0*(OE/endTem_475+1.0+OE1*endTem_475+OE2*endTem_475**2+OE3*endTem_475**3)
alphai_475=IA0*(IA/endTem_475+1.0+IA1*endTem_475+IA2*endTem_475**2+IA3*endTem_475**3)
alphao_475=OA0*(OA/endTem_475+1.0+OA1*endTem_475+OA2*endTem_475**2+OA3*endTem_475**3)
Ki_475=IK0*(IK/endTem_475+1.0+IK1*endTem_475+IK2*endTem_475**2+IK3*endTem_475**3)
Ko_475=OK0*(OK/endTem_475+1.0+OK1*endTem_475+OK2*endTem_475**2+OK3*endTem_475**3)
YS_475=MY0*(MY/endTem_475+1.0+MY1*endTem_475+MY2*endTem_475**2+MY3*endTem_475**3)
HS_475=MH0*(MH/endTem_475+1.0+MH1*endTem_475+MH2*endTem_475**2+MH3*endTem_475**3)

tanxine_475=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_475.append(Ei_475*Volumefenshu[i]+Eo_475*(1.0+(-Volumefenshu[i])))


pengzhang_475=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_475.append(alphai_475*Volumefenshu[i]+alphao_475*(1.0+(-Volumefenshu[i])))
    
    
daore_475=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_475.append(Ki_475*Volumefenshu[i]+Ko_475*(1.0+(-Volumefenshu[i])))
    

qufu_475_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_475_a.append(YS_475*(1.0+(-Volumefenshu[i])+Ei_475/Eo_475*Volumefenshu[i]))

   
qianghua_475=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_475.append(Ei_475*Volumefenshu[i]+HS_475*(1.0+(-Volumefenshu[i])))


qufu_475_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_475_b.append(qufu_475_a[i]+0.1*qianghua_475[i])


YS_475_a=YS_475   #Ti-6Al-4V
YS_475_b=YS_475_a+0.1*HS_475


#############�¶�500K�µĲ����������###############
Ei_500=IE0*(IE/endTem_500+1.0+IE1*endTem_500+IE2*endTem_500**2+IE3*endTem_500**3)
Eo_500=OE0*(OE/endTem_500+1.0+OE1*endTem_500+OE2*endTem_500**2+OE3*endTem_500**3)
alphai_500=IA0*(IA/endTem_500+1.0+IA1*endTem_500+IA2*endTem_500**2+IA3*endTem_500**3)
alphao_500=OA0*(OA/endTem_500+1.0+OA1*endTem_500+OA2*endTem_500**2+OA3*endTem_500**3)
Ki_500=IK0*(IK/endTem_500+1.0+IK1*endTem_500+IK2*endTem_500**2+IK3*endTem_500**3)
Ko_500=OK0*(OK/endTem_500+1.0+OK1*endTem_500+OK2*endTem_500**2+OK3*endTem_500**3)
YS_500=MY0*(MY/endTem_500+1.0+MY1*endTem_500+MY2*endTem_500**2+MY3*endTem_500**3)
HS_500=MH0*(MH/endTem_500+1.0+MH1*endTem_500+MH2*endTem_500**2+MH3*endTem_500**3)

tanxine_500=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_500.append(Ei_500*Volumefenshu[i]+Eo_500*(1.0+(-Volumefenshu[i])))


pengzhang_500=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_500.append(alphai_500*Volumefenshu[i]+alphao_500*(1.0+(-Volumefenshu[i])))
    

daore_500=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_500.append(Ki_500*Volumefenshu[i]+Ko_500*(1.0+(-Volumefenshu[i])))


qufu_500_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_500_a.append(YS_500*(1.0+(-Volumefenshu[i])+Ei_500/Eo_500*Volumefenshu[i]))

   
qianghua_500=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_500.append(Ei_500*Volumefenshu[i]+HS_500*(1.0+(-Volumefenshu[i])))


qufu_500_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_500_b.append(qufu_500_a[i]+0.1*qianghua_500[i])


YS_500_a=YS_500   #Ti-6Al-4V
YS_500_b=YS_500_a+0.1*HS_500


#############�¶�525K�µĲ����������###############
Ei_525=IE0*(IE/endTem_525+1.0+IE1*endTem_525+IE2*endTem_525**2+IE3*endTem_525**3)
Eo_525=OE0*(OE/endTem_525+1.0+OE1*endTem_525+OE2*endTem_525**2+OE3*endTem_525**3)
alphai_525=IA0*(IA/endTem_525+1.0+IA1*endTem_525+IA2*endTem_525**2+IA3*endTem_525**3)
alphao_525=OA0*(OA/endTem_525+1.0+OA1*endTem_525+OA2*endTem_525**2+OA3*endTem_525**3)
Ki_525=IK0*(IK/endTem_525+1.0+IK1*endTem_525+IK2*endTem_525**2+IK3*endTem_525**3)
Ko_525=OK0*(OK/endTem_525+1.0+OK1*endTem_525+OK2*endTem_525**2+OK3*endTem_525**3)
YS_525=MY0*(MY/endTem_525+1.0+MY1*endTem_525+MY2*endTem_525**2+MY3*endTem_525**3)
HS_525=MH0*(MH/endTem_525+1.0+MH1*endTem_525+MH2*endTem_525**2+MH3*endTem_525**3)

tanxine_525=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_525.append(Ei_525*Volumefenshu[i]+Eo_525*(1.0+(-Volumefenshu[i])))


pengzhang_525=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_525.append(alphai_525*Volumefenshu[i]+alphao_525*(1.0+(-Volumefenshu[i])))
    
    
daore_525=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_525.append(Ki_525*Volumefenshu[i]+Ko_525*(1.0+(-Volumefenshu[i])))


qufu_525_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_525_a.append(YS_525*(1.0+(-Volumefenshu[i])+Ei_525/Eo_525*Volumefenshu[i]))

   
qianghua_525=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_525.append(Ei_525*Volumefenshu[i]+HS_525*(1.0+(-Volumefenshu[i])))


qufu_525_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_525_b.append(qufu_525_a[i]+0.1*qianghua_525[i])
    

YS_525_a=YS_525   #Ti-6Al-4V
YS_525_b=YS_525_a+0.1*HS_525


#############�¶�550K�µĲ����������###############
Ei_550=IE0*(IE/endTem_550+1.0+IE1*endTem_550+IE2*endTem_550**2+IE3*endTem_550**3)
Eo_550=OE0*(OE/endTem_550+1.0+OE1*endTem_550+OE2*endTem_550**2+OE3*endTem_550**3)
alphai_550=IA0*(IA/endTem_550+1.0+IA1*endTem_550+IA2*endTem_550**2+IA3*endTem_550**3)
alphao_550=OA0*(OA/endTem_550+1.0+OA1*endTem_550+OA2*endTem_550**2+OA3*endTem_550**3)
Ki_550=IK0*(IK/endTem_550+1.0+IK1*endTem_550+IK2*endTem_550**2+IK3*endTem_550**3)
Ko_550=OK0*(OK/endTem_550+1.0+OK1*endTem_550+OK2*endTem_550**2+OK3*endTem_550**3)
YS_550=MY0*(MY/endTem_550+1.0+MY1*endTem_550+MY2*endTem_550**2+MY3*endTem_550**3)
HS_550=MH0*(MH/endTem_550+1.0+MH1*endTem_550+MH2*endTem_550**2+MH3*endTem_550**3)

tanxine_550=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_550.append(Ei_550*Volumefenshu[i]+Eo_550*(1.0+(-Volumefenshu[i])))


pengzhang_550=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_550.append(alphai_550*Volumefenshu[i]+alphao_550*(1.0+(-Volumefenshu[i])))
    
    
daore_550=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_550.append(Ki_550*Volumefenshu[i]+Ko_550*(1.0+(-Volumefenshu[i])))
    

qufu_550_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_550_a.append(YS_550*(1.0+(-Volumefenshu[i])+Ei_550/Eo_550*Volumefenshu[i]))

   
qianghua_550=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_550.append(Ei_550*Volumefenshu[i]+HS_550*(1.0+(-Volumefenshu[i])))


qufu_550_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_550_b.append(qufu_550_a[i]+0.1*qianghua_550[i])


YS_550_a=YS_550   #Ti-6Al-4V
YS_550_b=YS_550_a+0.1*HS_550


#############�¶�575K�µĲ����������###############
Ei_575=IE0*(IE/endTem_575+1.0+IE1*endTem_575+IE2*endTem_575**2+IE3*endTem_575**3)
Eo_575=OE0*(OE/endTem_575+1.0+OE1*endTem_575+OE2*endTem_575**2+OE3*endTem_575**3)
alphai_575=IA0*(IA/endTem_575+1.0+IA1*endTem_575+IA2*endTem_575**2+IA3*endTem_575**3)
alphao_575=OA0*(OA/endTem_575+1.0+OA1*endTem_575+OA2*endTem_575**2+OA3*endTem_575**3)
Ki_575=IK0*(IK/endTem_575+1.0+IK1*endTem_575+IK2*endTem_575**2+IK3*endTem_575**3)
Ko_575=OK0*(OK/endTem_575+1.0+OK1*endTem_575+OK2*endTem_575**2+OK3*endTem_575**3)
YS_575=MY0*(MY/endTem_575+1.0+MY1*endTem_575+MY2*endTem_575**2+MY3*endTem_575**3)
HS_575=MH0*(MH/endTem_575+1.0+MH1*endTem_575+MH2*endTem_575**2+MH3*endTem_575**3)

tanxine_575=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_575.append(Ei_575*Volumefenshu[i]+Eo_575*(1.0+(-Volumefenshu[i])))


pengzhang_575=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_575.append(alphai_575*Volumefenshu[i]+alphao_575*(1.0+(-Volumefenshu[i])))
    
    
daore_575=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_575.append(Ki_575*Volumefenshu[i]+Ko_575*(1.0+(-Volumefenshu[i])))


qufu_575_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_575_a.append(YS_575*(1.0+(-Volumefenshu[i])+Ei_575/Eo_575*Volumefenshu[i]))

   
qianghua_575=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_575.append(Ei_575*Volumefenshu[i]+HS_575*(1.0+(-Volumefenshu[i])))


qufu_575_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_575_b.append(qufu_575_a[i]+0.1*qianghua_575[i])


YS_575_a=YS_575   #Ti-6Al-4V
YS_575_b=YS_575_a+0.1*HS_575


#############�¶�600K�µĲ����������###############
Ei_600=IE0*(IE/endTem_600+1.0+IE1*endTem_600+IE2*endTem_600**2+IE3*endTem_600**3)
Eo_600=OE0*(OE/endTem_600+1.0+OE1*endTem_600+OE2*endTem_600**2+OE3*endTem_600**3)
alphai_600=IA0*(IA/endTem_600+1.0+IA1*endTem_600+IA2*endTem_600**2+IA3*endTem_600**3)
alphao_600=OA0*(OA/endTem_600+1.0+OA1*endTem_600+OA2*endTem_600**2+OA3*endTem_600**3)
Ki_600=IK0*(IK/endTem_600+1.0+IK1*endTem_600+IK2*endTem_600**2+IK3*endTem_600**3)
Ko_600=OK0*(OK/endTem_600+1.0+OK1*endTem_600+OK2*endTem_600**2+OK3*endTem_600**3)
YS_600=MY0*(MY/endTem_600+1.0+MY1*endTem_600+MY2*endTem_600**2+MY3*endTem_600**3)
HS_600=MH0*(MH/endTem_600+1.0+MH1*endTem_600+MH2*endTem_600**2+MH3*endTem_600**3)

tanxine_600=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_600.append(Ei_600*Volumefenshu[i]+Eo_600*(1.0+(-Volumefenshu[i])))


pengzhang_600=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_600.append(alphai_600*Volumefenshu[i]+alphao_600*(1.0+(-Volumefenshu[i])))
    
    
daore_600=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_600.append(Ki_600*Volumefenshu[i]+Ko_600*(1.0+(-Volumefenshu[i])))
    
    
qufu_600_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_600_a.append(YS_600*(1.0+(-Volumefenshu[i])+Ei_600/Eo_600*Volumefenshu[i]))

   
qianghua_600=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_600.append(Ei_600*Volumefenshu[i]+HS_600*(1.0+(-Volumefenshu[i])))


qufu_600_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_600_b.append(qufu_600_a[i]+0.1*qianghua_600[i])


YS_600_a=YS_600   #Ti-6Al-4V
YS_600_b=YS_600_a+0.1*HS_600


#############�¶�625K�µĲ����������###############
Ei_625=IE0*(IE/endTem_625+1.0+IE1*endTem_625+IE2*endTem_625**2+IE3*endTem_625**3)
Eo_625=OE0*(OE/endTem_625+1.0+OE1*endTem_625+OE2*endTem_625**2+OE3*endTem_625**3)
alphai_625=IA0*(IA/endTem_625+1.0+IA1*endTem_625+IA2*endTem_625**2+IA3*endTem_625**3)
alphao_625=OA0*(OA/endTem_625+1.0+OA1*endTem_625+OA2*endTem_625**2+OA3*endTem_625**3)
Ki_625=IK0*(IK/endTem_625+1.0+IK1*endTem_625+IK2*endTem_625**2+IK3*endTem_625**3)
Ko_625=OK0*(OK/endTem_625+1.0+OK1*endTem_625+OK2*endTem_625**2+OK3*endTem_625**3)
YS_625=MY0*(MY/endTem_625+1.0+MY1*endTem_625+MY2*endTem_625**2+MY3*endTem_625**3)
HS_625=MH0*(MH/endTem_625+1.0+MH1*endTem_625+MH2*endTem_625**2+MH3*endTem_625**3)

tanxine_625=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_625.append(Ei_625*Volumefenshu[i]+Eo_625*(1.0+(-Volumefenshu[i])))


pengzhang_625=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_625.append(alphai_625*Volumefenshu[i]+alphao_625*(1.0+(-Volumefenshu[i])))
    
    
daore_625=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_625.append(Ki_625*Volumefenshu[i]+Ko_625*(1.0+(-Volumefenshu[i])))
    
    
qufu_625_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_625_a.append(YS_625*(1.0+(-Volumefenshu[i])+Ei_625/Eo_625*Volumefenshu[i]))

   
qianghua_625=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_625.append(Ei_625*Volumefenshu[i]+HS_625*(1.0+(-Volumefenshu[i])))


qufu_625_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_625_b.append(qufu_625_a[i]+0.1*qianghua_625[i])
    
    
YS_625_a=YS_625   #Ti-6Al-4V
YS_625_b=YS_625_a+0.1*HS_625


#############�¶�650K�µĲ����������###############
Ei_650=IE0*(IE/endTem_650+1.0+IE1*endTem_650+IE2*endTem_650**2+IE3*endTem_650**3)
Eo_650=OE0*(OE/endTem_650+1.0+OE1*endTem_650+OE2*endTem_650**2+OE3*endTem_650**3)
alphai_650=IA0*(IA/endTem_650+1.0+IA1*endTem_650+IA2*endTem_650**2+IA3*endTem_650**3)
alphao_650=OA0*(OA/endTem_650+1.0+OA1*endTem_650+OA2*endTem_650**2+OA3*endTem_650**3)
Ki_650=IK0*(IK/endTem_650+1.0+IK1*endTem_650+IK2*endTem_650**2+IK3*endTem_650**3)
Ko_650=OK0*(OK/endTem_650+1.0+OK1*endTem_650+OK2*endTem_650**2+OK3*endTem_650**3)
YS_650=MY0*(MY/endTem_650+1.0+MY1*endTem_650+MY2*endTem_650**2+MY3*endTem_650**3)
HS_650=MH0*(MH/endTem_650+1.0+MH1*endTem_650+MH2*endTem_650**2+MH3*endTem_650**3)

tanxine_650=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_650.append(Ei_650*Volumefenshu[i]+Eo_650*(1.0+(-Volumefenshu[i])))


pengzhang_650=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_650.append(alphai_650*Volumefenshu[i]+alphao_650*(1.0+(-Volumefenshu[i])))
    
    
daore_650=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_650.append(Ki_650*Volumefenshu[i]+Ko_650*(1.0+(-Volumefenshu[i])))
    
    
qufu_650_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_650_a.append(YS_650*(1.0+(-Volumefenshu[i])+Ei_650/Eo_650*Volumefenshu[i]))

   
qianghua_650=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_650.append(Ei_650*Volumefenshu[i]+HS_650*(1.0+(-Volumefenshu[i])))


qufu_650_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_650_b.append(qufu_650_a[i]+0.1*qianghua_650[i])
    
    
YS_650_a=YS_650   #Ti-6Al-4V
YS_650_b=YS_650_a+0.1*HS_650

   
#############�¶�675K�µĲ����������###############
Ei_675=IE0*(IE/endTem_675+1.0+IE1*endTem_675+IE2*endTem_675**2+IE3*endTem_675**3)
Eo_675=OE0*(OE/endTem_675+1.0+OE1*endTem_675+OE2*endTem_675**2+OE3*endTem_675**3)
alphai_675=IA0*(IA/endTem_675+1.0+IA1*endTem_675+IA2*endTem_675**2+IA3*endTem_675**3)
alphao_675=OA0*(OA/endTem_675+1.0+OA1*endTem_675+OA2*endTem_675**2+OA3*endTem_675**3)
Ki_675=IK0*(IK/endTem_675+1.0+IK1*endTem_675+IK2*endTem_675**2+IK3*endTem_675**3)
Ko_675=OK0*(OK/endTem_675+1.0+OK1*endTem_675+OK2*endTem_675**2+OK3*endTem_675**3)
YS_675=MY0*(MY/endTem_675+1.0+MY1*endTem_675+MY2*endTem_675**2+MY3*endTem_675**3)
HS_675=MH0*(MH/endTem_675+1.0+MH1*endTem_675+MH2*endTem_675**2+MH3*endTem_675**3)

tanxine_675=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_675.append(Ei_675*Volumefenshu[i]+Eo_675*(1.0+(-Volumefenshu[i])))


pengzhang_675=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_675.append(alphai_675*Volumefenshu[i]+alphao_675*(1.0+(-Volumefenshu[i])))
    
    
daore_675=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_675.append(Ki_675*Volumefenshu[i]+Ko_675*(1.0+(-Volumefenshu[i])))
    
    
qufu_675_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_675_a.append(YS_675*(1.0+(-Volumefenshu[i])+Ei_675/Eo_675*Volumefenshu[i]))

   
qianghua_675=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_675.append(Ei_675*Volumefenshu[i]+HS_675*(1.0+(-Volumefenshu[i])))


qufu_675_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_675_b.append(qufu_675_a[i]+0.1*qianghua_675[i])
    
    
YS_675_a=YS_675   #Ti-6Al-4V
YS_675_b=YS_675_a+0.1*HS_675


#############�¶�700K�µĲ����������###############
Ei_700=IE0*(IE/endTem_700+1.0+IE1*endTem_700+IE2*endTem_700**2+IE3*endTem_700**3)
Eo_700=OE0*(OE/endTem_700+1.0+OE1*endTem_700+OE2*endTem_700**2+OE3*endTem_700**3)
alphai_700=IA0*(IA/endTem_700+1.0+IA1*endTem_700+IA2*endTem_700**2+IA3*endTem_700**3)
alphao_700=OA0*(OA/endTem_700+1.0+OA1*endTem_700+OA2*endTem_700**2+OA3*endTem_700**3)
Ki_700=IK0*(IK/endTem_700+1.0+IK1*endTem_700+IK2*endTem_700**2+IK3*endTem_700**3)
Ko_700=OK0*(OK/endTem_700+1.0+OK1*endTem_700+OK2*endTem_700**2+OK3*endTem_700**3)
YS_700=MY0*(MY/endTem_700+1.0+MY1*endTem_700+MY2*endTem_700**2+MY3*endTem_700**3)
HS_700=MH0*(MH/endTem_700+1.0+MH1*endTem_700+MH2*endTem_700**2+MH3*endTem_700**3)

tanxine_700=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_700.append(Ei_700*Volumefenshu[i]+Eo_700*(1.0+(-Volumefenshu[i])))


pengzhang_700=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_700.append(alphai_700*Volumefenshu[i]+alphao_700*(1.0+(-Volumefenshu[i])))
    
    
daore_700=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_700.append(Ki_700*Volumefenshu[i]+Ko_700*(1.0+(-Volumefenshu[i])))
    
    
qufu_700_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_700_a.append(YS_700*(1.0+(-Volumefenshu[i])+Ei_700/Eo_700*Volumefenshu[i]))

   
qianghua_700=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_700.append(Ei_700*Volumefenshu[i]+HS_700*(1.0+(-Volumefenshu[i])))


qufu_700_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_700_b.append(qufu_700_a[i]+0.1*qianghua_700[i])
    
    
YS_700_a=YS_700   #Ti-6Al-4V
YS_700_b=YS_700_a+0.1*HS_700


#############�¶�725K�µĲ����������###############
Ei_725=IE0*(IE/endTem_725+1.0+IE1*endTem_725+IE2*endTem_725**2+IE3*endTem_725**3)
Eo_725=OE0*(OE/endTem_725+1.0+OE1*endTem_725+OE2*endTem_725**2+OE3*endTem_725**3)
alphai_725=IA0*(IA/endTem_725+1.0+IA1*endTem_725+IA2*endTem_725**2+IA3*endTem_725**3)
alphao_725=OA0*(OA/endTem_725+1.0+OA1*endTem_725+OA2*endTem_725**2+OA3*endTem_725**3)
Ki_725=IK0*(IK/endTem_725+1.0+IK1*endTem_725+IK2*endTem_725**2+IK3*endTem_725**3)
Ko_725=OK0*(OK/endTem_725+1.0+OK1*endTem_725+OK2*endTem_725**2+OK3*endTem_725**3)
YS_725=MY0*(MY/endTem_725+1.0+MY1*endTem_725+MY2*endTem_725**2+MY3*endTem_725**3)
HS_725=MH0*(MH/endTem_725+1.0+MH1*endTem_725+MH2*endTem_725**2+MH3*endTem_725**3)

tanxine_725=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_725.append(Ei_725*Volumefenshu[i]+Eo_725*(1.0+(-Volumefenshu[i])))


pengzhang_725=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_725.append(alphai_725*Volumefenshu[i]+alphao_725*(1.0+(-Volumefenshu[i])))
    
    
daore_725=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_725.append(Ki_725*Volumefenshu[i]+Ko_725*(1.0+(-Volumefenshu[i])))
    
    
qufu_725_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_725_a.append(YS_725*(1.0+(-Volumefenshu[i])+Ei_725/Eo_725*Volumefenshu[i]))

   
qianghua_725=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_725.append(Ei_725*Volumefenshu[i]+HS_725*(1.0+(-Volumefenshu[i])))


qufu_725_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_725_b.append(qufu_725_a[i]+0.1*qianghua_725[i])
    
    
YS_725_a=YS_725   #Ti-6Al-4V
YS_725_b=YS_725_a+0.1*HS_725


#############�¶�750K�µĲ����������###############
Ei_750=IE0*(IE/endTem_750+1.0+IE1*endTem_750+IE2*endTem_750**2+IE3*endTem_750**3)
Eo_750=OE0*(OE/endTem_750+1.0+OE1*endTem_750+OE2*endTem_750**2+OE3*endTem_750**3)
alphai_750=IA0*(IA/endTem_750+1.0+IA1*endTem_750+IA2*endTem_750**2+IA3*endTem_750**3)
alphao_750=OA0*(OA/endTem_750+1.0+OA1*endTem_750+OA2*endTem_750**2+OA3*endTem_750**3)
Ki_750=IK0*(IK/endTem_750+1.0+IK1*endTem_750+IK2*endTem_750**2+IK3*endTem_750**3)
Ko_750=OK0*(OK/endTem_750+1.0+OK1*endTem_750+OK2*endTem_750**2+OK3*endTem_750**3)
YS_750=MY0*(MY/endTem_750+1.0+MY1*endTem_750+MY2*endTem_750**2+MY3*endTem_750**3)
HS_750=MH0*(MH/endTem_750+1.0+MH1*endTem_750+MH2*endTem_750**2+MH3*endTem_750**3)

tanxine_750=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_750.append(Ei_750*Volumefenshu[i]+Eo_750*(1.0+(-Volumefenshu[i])))


pengzhang_750=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_750.append(alphai_750*Volumefenshu[i]+alphao_750*(1.0+(-Volumefenshu[i])))
    
    
daore_750=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_750.append(Ki_750*Volumefenshu[i]+Ko_750*(1.0+(-Volumefenshu[i])))
    
    
qufu_750_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_750_a.append(YS_750*(1.0+(-Volumefenshu[i])+Ei_750/Eo_750*Volumefenshu[i]))

   
qianghua_750=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_750.append(Ei_750*Volumefenshu[i]+HS_750*(1.0+(-Volumefenshu[i])))


qufu_750_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_750_b.append(qufu_750_a[i]+0.1*qianghua_750[i])
    
    
YS_750_a=YS_750   #Ti-6Al-4V
YS_750_b=YS_750_a+0.1*HS_750


#############�¶�775K�µĲ����������###############
Ei_775=IE0*(IE/endTem_775+1.0+IE1*endTem_775+IE2*endTem_775**2+IE3*endTem_775**3)
Eo_775=OE0*(OE/endTem_775+1.0+OE1*endTem_775+OE2*endTem_775**2+OE3*endTem_775**3)
alphai_775=IA0*(IA/endTem_775+1.0+IA1*endTem_775+IA2*endTem_775**2+IA3*endTem_775**3)
alphao_775=OA0*(OA/endTem_775+1.0+OA1*endTem_775+OA2*endTem_775**2+OA3*endTem_775**3)
Ki_775=IK0*(IK/endTem_775+1.0+IK1*endTem_775+IK2*endTem_775**2+IK3*endTem_775**3)
Ko_775=OK0*(OK/endTem_775+1.0+OK1*endTem_775+OK2*endTem_775**2+OK3*endTem_775**3)
YS_775=MY0*(MY/endTem_775+1.0+MY1*endTem_775+MY2*endTem_775**2+MY3*endTem_775**3)
HS_775=MH0*(MH/endTem_775+1.0+MH1*endTem_775+MH2*endTem_775**2+MH3*endTem_775**3)

tanxine_775=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_775.append(Ei_775*Volumefenshu[i]+Eo_775*(1.0+(-Volumefenshu[i])))


pengzhang_775=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_775.append(alphai_775*Volumefenshu[i]+alphao_775*(1.0+(-Volumefenshu[i])))
    
    
daore_775=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_775.append(Ki_775*Volumefenshu[i]+Ko_775*(1.0+(-Volumefenshu[i])))
    
    
qufu_775_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_775_a.append(YS_775*(1.0+(-Volumefenshu[i])+Ei_775/Eo_775*Volumefenshu[i]))

   
qianghua_775=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_775.append(Ei_775*Volumefenshu[i]+HS_775*(1.0+(-Volumefenshu[i])))


qufu_775_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_775_b.append(qufu_775_a[i]+0.1*qianghua_775[i])
    
    
YS_775_a=YS_775   #Ti-6Al-4V
YS_775_b=YS_775_a+0.1*HS_775


#############�¶�800K�µĲ����������###############
Ei_800=IE0*(IE/endTem_800+1.0+IE1*endTem_800+IE2*endTem_800**2+IE3*endTem_800**3)
Eo_800=OE0*(OE/endTem_800+1.0+OE1*endTem_800+OE2*endTem_800**2+OE3*endTem_800**3)
alphai_800=IA0*(IA/endTem_800+1.0+IA1*endTem_800+IA2*endTem_800**2+IA3*endTem_800**3)
alphao_800=OA0*(OA/endTem_800+1.0+OA1*endTem_800+OA2*endTem_800**2+OA3*endTem_800**3)
Ki_800=IK0*(IK/endTem_800+1.0+IK1*endTem_800+IK2*endTem_800**2+IK3*endTem_800**3)
Ko_800=OK0*(OK/endTem_800+1.0+OK1*endTem_800+OK2*endTem_800**2+OK3*endTem_800**3)
YS_800=MY0*(MY/endTem_800+1.0+MY1*endTem_800+MY2*endTem_800**2+MY3*endTem_800**3)
HS_800=MH0*(MH/endTem_800+1.0+MH1*endTem_800+MH2*endTem_800**2+MH3*endTem_800**3)

tanxine_800=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_800.append(Ei_800*Volumefenshu[i]+Eo_800*(1.0+(-Volumefenshu[i])))


pengzhang_800=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_800.append(alphai_800*Volumefenshu[i]+alphao_800*(1.0+(-Volumefenshu[i])))
    
    
daore_800=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_800.append(Ki_800*Volumefenshu[i]+Ko_800*(1.0+(-Volumefenshu[i])))
    
    
qufu_800_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_800_a.append(YS_800*(1.0+(-Volumefenshu[i])+Ei_800/Eo_800*Volumefenshu[i]))

   
qianghua_800=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_800.append(Ei_800*Volumefenshu[i]+HS_800*(1.0+(-Volumefenshu[i])))


qufu_800_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_800_b.append(qufu_800_a[i]+0.1*qianghua_800[i])
    
    
YS_800_a=YS_800   #Ti-6Al-4V
YS_800_b=YS_800_a+0.1*HS_800


#############�¶�825K�µĲ����������###############
Ei_825=IE0*(IE/endTem_825+1.0+IE1*endTem_825+IE2*endTem_825**2+IE3*endTem_825**3)
Eo_825=OE0*(OE/endTem_825+1.0+OE1*endTem_825+OE2*endTem_825**2+OE3*endTem_825**3)
alphai_825=IA0*(IA/endTem_825+1.0+IA1*endTem_825+IA2*endTem_825**2+IA3*endTem_825**3)
alphao_825=OA0*(OA/endTem_825+1.0+OA1*endTem_825+OA2*endTem_825**2+OA3*endTem_825**3)
Ki_825=IK0*(IK/endTem_825+1.0+IK1*endTem_825+IK2*endTem_825**2+IK3*endTem_825**3)
Ko_825=OK0*(OK/endTem_825+1.0+OK1*endTem_825+OK2*endTem_825**2+OK3*endTem_825**3)
YS_825=MY0*(MY/endTem_825+1.0+MY1*endTem_825+MY2*endTem_825**2+MY3*endTem_825**3)
HS_825=MH0*(MH/endTem_825+1.0+MH1*endTem_825+MH2*endTem_825**2+MH3*endTem_825**3)

tanxine_825=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_825.append(Ei_825*Volumefenshu[i]+Eo_825*(1.0+(-Volumefenshu[i])))


pengzhang_825=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_825.append(alphai_825*Volumefenshu[i]+alphao_825*(1.0+(-Volumefenshu[i])))
    
    
daore_825=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_825.append(Ki_825*Volumefenshu[i]+Ko_825*(1.0+(-Volumefenshu[i])))
    
    
qufu_825_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_825_a.append(YS_825*(1.0+(-Volumefenshu[i])+Ei_825/Eo_825*Volumefenshu[i]))

   
qianghua_825=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_825.append(Ei_825*Volumefenshu[i]+HS_825*(1.0+(-Volumefenshu[i])))


qufu_825_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_825_b.append(qufu_825_a[i]+0.1*qianghua_825[i])
    
    
YS_825_a=YS_825   #Ti-6Al-4V
YS_825_b=YS_825_a+0.1*HS_825


#############�¶�850K�µĲ����������###############
Ei_850=IE0*(IE/endTem_850+1.0+IE1*endTem_850+IE2*endTem_850**2+IE3*endTem_850**3)
Eo_850=OE0*(OE/endTem_850+1.0+OE1*endTem_850+OE2*endTem_850**2+OE3*endTem_850**3)
alphai_850=IA0*(IA/endTem_850+1.0+IA1*endTem_850+IA2*endTem_850**2+IA3*endTem_850**3)
alphao_850=OA0*(OA/endTem_850+1.0+OA1*endTem_850+OA2*endTem_850**2+OA3*endTem_850**3)
Ki_850=IK0*(IK/endTem_850+1.0+IK1*endTem_850+IK2*endTem_850**2+IK3*endTem_850**3)
Ko_850=OK0*(OK/endTem_850+1.0+OK1*endTem_850+OK2*endTem_850**2+OK3*endTem_850**3)
YS_850=MY0*(MY/endTem_850+1.0+MY1*endTem_850+MY2*endTem_850**2+MY3*endTem_850**3)
HS_850=MH0*(MH/endTem_850+1.0+MH1*endTem_850+MH2*endTem_850**2+MH3*endTem_850**3)

tanxine_850=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_850.append(Ei_850*Volumefenshu[i]+Eo_850*(1.0+(-Volumefenshu[i])))


pengzhang_850=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_850.append(alphai_850*Volumefenshu[i]+alphao_850*(1.0+(-Volumefenshu[i])))
    
    
daore_850=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_850.append(Ki_850*Volumefenshu[i]+Ko_850*(1.0+(-Volumefenshu[i])))
    
    
qufu_850_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_850_a.append(YS_850*(1.0+(-Volumefenshu[i])+Ei_850/Eo_850*Volumefenshu[i]))

   
qianghua_850=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_850.append(Ei_850*Volumefenshu[i]+HS_850*(1.0+(-Volumefenshu[i])))


qufu_850_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_850_b.append(qufu_850_a[i]+0.1*qianghua_850[i])
    
    
YS_850_a=YS_850   #Ti-6Al-4V
YS_850_b=YS_850_a+0.1*HS_850


#############�¶�875K�µĲ����������###############
Ei_875=IE0*(IE/endTem_875+1.0+IE1*endTem_875+IE2*endTem_875**2+IE3*endTem_875**3)
Eo_875=OE0*(OE/endTem_875+1.0+OE1*endTem_875+OE2*endTem_875**2+OE3*endTem_875**3)
alphai_875=IA0*(IA/endTem_875+1.0+IA1*endTem_875+IA2*endTem_875**2+IA3*endTem_875**3)
alphao_875=OA0*(OA/endTem_875+1.0+OA1*endTem_875+OA2*endTem_875**2+OA3*endTem_875**3)
Ki_875=IK0*(IK/endTem_875+1.0+IK1*endTem_875+IK2*endTem_875**2+IK3*endTem_875**3)
Ko_875=OK0*(OK/endTem_875+1.0+OK1*endTem_875+OK2*endTem_875**2+OK3*endTem_875**3)
YS_875=MY0*(MY/endTem_875+1.0+MY1*endTem_875+MY2*endTem_875**2+MY3*endTem_875**3)
HS_875=MH0*(MH/endTem_875+1.0+MH1*endTem_875+MH2*endTem_875**2+MH3*endTem_875**3)

tanxine_875=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_875.append(Ei_875*Volumefenshu[i]+Eo_875*(1.0+(-Volumefenshu[i])))


pengzhang_875=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_875.append(alphai_875*Volumefenshu[i]+alphao_875*(1.0+(-Volumefenshu[i])))
    
    
daore_875=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_875.append(Ki_875*Volumefenshu[i]+Ko_875*(1.0+(-Volumefenshu[i])))
    
    
qufu_875_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_875_a.append(YS_875*(1.0+(-Volumefenshu[i])+Ei_875/Eo_875*Volumefenshu[i]))

   
qianghua_875=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_875.append(Ei_875*Volumefenshu[i]+HS_875*(1.0+(-Volumefenshu[i])))


qufu_875_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_875_b.append(qufu_875_a[i]+0.1*qianghua_875[i])
    
    
YS_875_a=YS_875   #Ti-6Al-4V
YS_875_b=YS_875_a+0.1*HS_875


#############�¶�900K�µĲ����������###############
Ei_900=IE0*(IE/endTem_900+1.0+IE1*endTem_900+IE2*endTem_900**2+IE3*endTem_900**3)
Eo_900=OE0*(OE/endTem_900+1.0+OE1*endTem_900+OE2*endTem_900**2+OE3*endTem_900**3)
alphai_900=IA0*(IA/endTem_900+1.0+IA1*endTem_900+IA2*endTem_900**2+IA3*endTem_900**3)
alphao_900=OA0*(OA/endTem_900+1.0+OA1*endTem_900+OA2*endTem_900**2+OA3*endTem_900**3)
Ki_900=IK0*(IK/endTem_900+1.0+IK1*endTem_900+IK2*endTem_900**2+IK3*endTem_900**3)
Ko_900=OK0*(OK/endTem_900+1.0+OK1*endTem_900+OK2*endTem_900**2+OK3*endTem_900**3)
YS_900=MY0*(MY/endTem_900+1.0+MY1*endTem_900+MY2*endTem_900**2+MY3*endTem_900**3)
HS_900=MH0*(MH/endTem_900+1.0+MH1*endTem_900+MH2*endTem_900**2+MH3*endTem_900**3)

tanxine_900=[]    #����ģ�������¶ȡ������йأ�N/mm^2
for i in range(0, cengshu): 
    tanxine_900.append(Ei_900*Volumefenshu[i]+Eo_900*(1.0+(-Volumefenshu[i])))


pengzhang_900=[]  #������ϵ�������¶ȡ������йأ�1/K
for i in range(0, cengshu):         
    pengzhang_900.append(alphai_900*Volumefenshu[i]+alphao_900*(1.0+(-Volumefenshu[i])))
    
    
daore_900=[]      #����ϵ��(�ȵ���)�����¶ȡ������йأ�W/m*K
for i in range(0, cengshu):         
    daore_900.append(Ki_900*Volumefenshu[i]+Ko_900*(1.0+(-Volumefenshu[i])))
    
    
qufu_900_a=[]     #��������(ʼ�㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_900_a.append(YS_900*(1.0+(-Volumefenshu[i])+Ei_900/Eo_900*Volumefenshu[i]))

   
qianghua_900=[]   #ǿ��ģ�������¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qianghua_900.append(Ei_900*Volumefenshu[i]+HS_900*(1.0+(-Volumefenshu[i])))


qufu_900_b=[]     #��������(�յ㣩�����¶ȡ������йأ�MPa
for i in range(0, cengshu):         
    qufu_900_b.append(qufu_900_a[i]+0.1*qianghua_900[i])
    
    
YS_900_a=YS_900   #Ti-6Al-4V
YS_900_b=YS_900_a+0.1*HS_900


#############################################################
###############################################################
###################################################################
Mdb()
mdb.Model(name=modelname,modelType=STANDARD_EXPLICIT)     #modelname='aa'

m=mdb.models[modelname]
session.journalOptions.setValues(replayGeometry=COORDINATE,
    recoverGeometry=COORDINATE)        #�����������ǿɸ�������ѡ�����
r=m.rootAssembly

#################
for i in range(1, cengshu+1): 
    matname='O-I-'+str(i)
    mdb.models[modelname].Material(name=matname)         #����-�մ�, O-I-1�� O-I-cengshu�����ڲ㵽���
    mdb.models[modelname].materials[matname].Density(table=((midu[i-1], ), ))
    mdb.models[modelname].materials[matname].Elastic(temperatureDependency=ON,
        table=((tanxine_300[i-1],posongb,endTem_300), (tanxine_325[i-1],posongb,endTem_325), (tanxine_350[i-1],posongb,endTem_350), (tanxine_375[i-1],posongb,endTem_375), 
        (tanxine_400[i-1],posongb,endTem_400), (tanxine_425[i-1],posongb,endTem_425), (tanxine_450[i-1],posongb,endTem_450), (tanxine_475[i-1],posongb,endTem_475), 
        (tanxine_500[i-1],posongb,endTem_500), (tanxine_525[i-1],posongb,endTem_525), (tanxine_550[i-1],posongb,endTem_550), (tanxine_575[i-1],posongb,endTem_575), 
        (tanxine_600[i-1],posongb,endTem_600), (tanxine_625[i-1],posongb,endTem_625), (tanxine_650[i-1],posongb,endTem_650), (tanxine_675[i-1],posongb,endTem_675), 
        (tanxine_700[i-1],posongb,endTem_700), (tanxine_725[i-1],posongb,endTem_725), (tanxine_750[i-1],posongb,endTem_750), (tanxine_775[i-1],posongb,endTem_775), 
        (tanxine_800[i-1],posongb,endTem_800), (tanxine_825[i-1],posongb,endTem_825), (tanxine_850[i-1],posongb,endTem_850), (tanxine_875[i-1],posongb,endTem_875), 
        (tanxine_900[i-1],posongb,endTem_900)))
    mdb.models[modelname].materials[matname].Plastic(temperatureDependency=ON,table=((qufu_300_a[i-1], 0.0, endTem_300), (qufu_300_b[i-1], 0.1, endTem_300), (qufu_325_a[i-1], 0.0, endTem_325), (qufu_325_b[i-1], 0.1, endTem_325), 
        (qufu_350_a[i-1], 0.0, endTem_350), (qufu_350_b[i-1], 0.1, endTem_350), (qufu_375_a[i-1], 0.0, endTem_375), (qufu_375_b[i-1], 0.1, endTem_375), (qufu_400_a[i-1], 0.0, endTem_400), (qufu_400_b[i-1], 0.1, endTem_400), 
        (qufu_425_a[i-1], 0.0, endTem_425), (qufu_425_b[i-1], 0.1, endTem_425), (qufu_450_a[i-1], 0.0, endTem_450), (qufu_450_b[i-1], 0.1, endTem_450), (qufu_475_a[i-1], 0.0, endTem_475), (qufu_475_b[i-1], 0.1, endTem_475), 
        (qufu_500_a[i-1], 0.0, endTem_500), (qufu_500_b[i-1], 0.1, endTem_500), (qufu_525_a[i-1], 0.0, endTem_525), (qufu_525_b[i-1], 0.1, endTem_525), (qufu_550_a[i-1], 0.0, endTem_550), (qufu_550_b[i-1], 0.1, endTem_550), 
        (qufu_575_a[i-1], 0.0, endTem_575), (qufu_575_b[i-1], 0.1, endTem_575), (qufu_600_a[i-1], 0.0, endTem_600), (qufu_600_b[i-1], 0.1, endTem_600), (qufu_625_a[i-1], 0.0, endTem_625), (qufu_625_b[i-1], 0.1, endTem_625), 
        (qufu_650_a[i-1], 0.0, endTem_650), (qufu_650_b[i-1], 0.1, endTem_650), (qufu_675_a[i-1], 0.0, endTem_675), (qufu_675_b[i-1], 0.1, endTem_675), (qufu_700_a[i-1], 0.0, endTem_700), (qufu_700_b[i-1], 0.1, endTem_700), 
        (qufu_725_a[i-1], 0.0, endTem_725), (qufu_725_b[i-1], 0.1, endTem_725), (qufu_750_a[i-1], 0.0, endTem_750), (qufu_750_b[i-1], 0.1, endTem_750), (qufu_775_a[i-1], 0.0, endTem_775), (qufu_775_b[i-1], 0.1, endTem_775), 
        (qufu_800_a[i-1], 0.0, endTem_800), (qufu_800_b[i-1], 0.1, endTem_800), (qufu_825_a[i-1], 0.0, endTem_825), (qufu_825_b[i-1], 0.1, endTem_825), (qufu_850_a[i-1], 0.0, endTem_850), (qufu_850_b[i-1], 0.1, endTem_850), 
        (qufu_875_a[i-1], 0.0, endTem_875), (qufu_875_b[i-1], 0.1, endTem_875), (qufu_900_a[i-1], 0.0, endTem_900), (qufu_900_b[i-1], 0.1, endTem_900), ))
    mdb.models[modelname].materials[matname].Expansion(temperatureDependency=ON, 
        table=((pengzhang_300[i-1], endTem_300), (pengzhang_325[i-1], endTem_325), (pengzhang_350[i-1], endTem_350), (pengzhang_375[i-1], endTem_375), 
        (pengzhang_400[i-1], endTem_400), (pengzhang_425[i-1], endTem_425), (pengzhang_450[i-1], endTem_450), (pengzhang_475[i-1], endTem_475), 
        (pengzhang_500[i-1], endTem_500), (pengzhang_525[i-1], endTem_525), (pengzhang_550[i-1], endTem_550), (pengzhang_575[i-1], endTem_575), 
        (pengzhang_600[i-1], endTem_600), (pengzhang_625[i-1], endTem_625), (pengzhang_650[i-1], endTem_650), (pengzhang_675[i-1], endTem_675), 
        (pengzhang_700[i-1], endTem_700), (pengzhang_725[i-1], endTem_725), (pengzhang_750[i-1], endTem_750), (pengzhang_775[i-1], endTem_775), 
        (pengzhang_800[i-1], endTem_800), (pengzhang_825[i-1], endTem_825), (pengzhang_850[i-1], endTem_850), (pengzhang_875[i-1], endTem_875), 
        (pengzhang_900[i-1], endTem_900)), zero=298.15)
    mdb.models[modelname].materials[matname].Conductivity(temperatureDependency=ON,
        table=((daore_300[i-1], endTem_300), (daore_325[i-1], endTem_325), (daore_350[i-1], endTem_350), (daore_375[i-1], endTem_375), 
        (daore_400[i-1], endTem_400), (daore_425[i-1], endTem_425), (daore_450[i-1], endTem_450), (daore_475[i-1], endTem_475), 
        (daore_500[i-1], endTem_500), (daore_525[i-1], endTem_525), (daore_550[i-1], endTem_550), (daore_575[i-1], endTem_575), 
        (daore_600[i-1], endTem_600), (daore_625[i-1], endTem_625), (daore_650[i-1], endTem_650), (daore_675[i-1], endTem_675), 
        (daore_700[i-1], endTem_700), (daore_725[i-1], endTem_725), (daore_750[i-1], endTem_750), (daore_775[i-1], endTem_775), 
        (daore_800[i-1], endTem_800), (daore_825[i-1], endTem_825), (daore_850[i-1], endTem_850), (daore_875[i-1], endTem_875), 
        (daore_900[i-1], endTem_900)))



mdb.models[modelname].Material(name='AlO')             #�մɣ�Al2O3
mdb.models[modelname].materials['AlO'].Density(table=((Omidu, ), ))
mdb.models[modelname].materials['AlO'].Elastic(temperatureDependency=ON,
    table=((Eo_300,Oposongb,endTem_300), (Eo_325,Oposongb,endTem_325), (Eo_350,Oposongb,endTem_350), (Eo_375,Oposongb,endTem_375), 
    (Eo_400,Oposongb,endTem_400), (Eo_425,Oposongb,endTem_425), (Eo_450,Oposongb,endTem_450), (Eo_475,Oposongb,endTem_475), 
    (Eo_500,Oposongb,endTem_500), (Eo_525,Oposongb,endTem_525), (Eo_550,Oposongb,endTem_550), (Eo_575,Oposongb,endTem_575), 
    (Eo_600,Oposongb,endTem_600), (Eo_625,Oposongb,endTem_625), (Eo_650,Oposongb,endTem_650), (Eo_675,Oposongb,endTem_675), 
    (Eo_700,Oposongb,endTem_700), (Eo_725,Oposongb,endTem_725), (Eo_750,Oposongb,endTem_750), (Eo_775,Oposongb,endTem_775), 
    (Eo_800,Oposongb,endTem_800), (Eo_825,Oposongb,endTem_825), (Eo_850,Oposongb,endTem_850), (Eo_875,Oposongb,endTem_875), 
    (Eo_900,Oposongb,endTem_900)))
#mdb.models[modelname].materials['AlO'].Plastic(table=((YS, 0.0), ))  #�մɲ���������
mdb.models[modelname].materials['AlO'].Expansion(temperatureDependency=ON, 
    table=((alphao_300, endTem_300), (alphao_325, endTem_325), (alphao_350, endTem_350), (alphao_375, endTem_375), 
    (alphao_400, endTem_400), (alphao_425, endTem_425), (alphao_450, endTem_450), (alphao_475, endTem_475), 
    (alphao_500, endTem_500), (alphao_525, endTem_525), (alphao_550, endTem_550), (alphao_575, endTem_575), 
    (alphao_600, endTem_600), (alphao_625, endTem_625), (alphao_650, endTem_650), (alphao_675, endTem_675), 
    (alphao_700, endTem_700), (alphao_725, endTem_725), (alphao_750, endTem_750), (alphao_775, endTem_775), 
    (alphao_800, endTem_800), (alphao_825, endTem_825), (alphao_850, endTem_850), (alphao_875, endTem_875), 
    (alphao_900, endTem_900)), zero=298.15)
mdb.models[modelname].materials['AlO'].Conductivity(temperatureDependency=ON,
     table=((Ko_300, endTem_300), (Ko_325, endTem_325), (Ko_350, endTem_350), (Ko_375, endTem_375), 
     (Ko_400, endTem_400), (Ko_425, endTem_425), (Ko_450, endTem_450), (Ko_475, endTem_475), 
     (Ko_500, endTem_500), (Ko_525, endTem_525), (Ko_550, endTem_550), (Ko_575, endTem_575), 
     (Ko_600, endTem_600), (Ko_625, endTem_625), (Ko_650, endTem_650), (Ko_675, endTem_675), 
     (Ko_700, endTem_700), (Ko_725, endTem_725), (Ko_750, endTem_750), (Ko_775, endTem_775), 
     (Ko_800, endTem_800), (Ko_825, endTem_825), (Ko_850, endTem_850), (Ko_875, endTem_875), 
     (Ko_900, endTem_900)))
     
mdb.models[modelname].Material(name='Ti')            #������Ti-6Al-4V
mdb.models[modelname].materials['Ti'].Density(table=((Imidu, ), ))
mdb.models[modelname].materials['Ti'].Elastic(temperatureDependency=ON,
    table=((Ei_300,Iposongb,endTem_300), (Ei_325,Iposongb,endTem_325), (Ei_350,Iposongb,endTem_350), (Ei_375,Iposongb,endTem_375), 
    (Ei_400,Iposongb,endTem_400), (Ei_425,Iposongb,endTem_425), (Ei_450,Iposongb,endTem_450), (Ei_475,Iposongb,endTem_475), 
    (Ei_500,Iposongb,endTem_500), (Ei_525,Iposongb,endTem_525), (Ei_550,Iposongb,endTem_550), (Ei_575,Iposongb,endTem_575), 
    (Ei_600,Iposongb,endTem_600), (Ei_625,Iposongb,endTem_625), (Ei_650,Iposongb,endTem_650), (Ei_675,Iposongb,endTem_675), 
    (Ei_700,Iposongb,endTem_700), (Ei_725,Iposongb,endTem_725), (Ei_750,Iposongb,endTem_750), (Ei_775,Iposongb,endTem_775), 
    (Ei_800,Iposongb,endTem_800), (Ei_825,Iposongb,endTem_825), (Ei_850,Iposongb,endTem_850), (Ei_875,Iposongb,endTem_875), 
    (Ei_900,Iposongb,endTem_900)))
mdb.models[modelname].materials['Ti'].Plastic(temperatureDependency=ON,table=((YS_300_a, 0.0, endTem_300), (YS_300_b, 0.1, endTem_300), (YS_325_a, 0.0, endTem_325), (YS_325_b, 0.1, endTem_325), 
    (YS_350_a, 0.0, endTem_350), (YS_350_b, 0.1, endTem_350), (YS_375_a, 0.0, endTem_375), (YS_375_b, 0.1, endTem_375), (YS_400_a, 0.0, endTem_400), (YS_400_b, 0.1, endTem_400), 
    (YS_425_a, 0.0, endTem_425), (YS_425_b, 0.1, endTem_425), (YS_450_a, 0.0, endTem_450), (YS_450_b, 0.1, endTem_450), (YS_475_a, 0.0, endTem_475), (YS_475_b, 0.1, endTem_475), 
    (YS_500_a, 0.0, endTem_500), (YS_500_b, 0.1, endTem_500), (YS_525_a, 0.0, endTem_525), (YS_525_b, 0.1, endTem_525), (YS_550_a, 0.0, endTem_550), (YS_550_b, 0.1, endTem_550), 
    (YS_575_a, 0.0, endTem_575), (YS_575_b, 0.1, endTem_575), (YS_600_a, 0.0, endTem_600), (YS_600_b, 0.1, endTem_600), (YS_625_a, 0.0, endTem_625), (YS_625_b, 0.1, endTem_625), 
    (YS_650_a, 0.0, endTem_650), (YS_650_b, 0.1, endTem_650), (YS_675_a, 0.0, endTem_675), (YS_675_b, 0.1, endTem_675), (YS_700_a, 0.0, endTem_700), (YS_700_b, 0.1, endTem_700), 
    (YS_725_a, 0.0, endTem_725), (YS_725_b, 0.1, endTem_725), (YS_750_a, 0.0, endTem_750), (YS_750_b, 0.1, endTem_750), (YS_775_a, 0.0, endTem_775), (YS_775_b, 0.1, endTem_775), 
    (YS_800_a, 0.0, endTem_800), (YS_800_b, 0.1, endTem_800), (YS_825_a, 0.0, endTem_825), (YS_825_b, 0.1, endTem_825), (YS_850_a, 0.0, endTem_850), (YS_850_b, 0.1, endTem_850), 
    (YS_875_a, 0.0, endTem_875), (YS_875_b, 0.1, endTem_875), (YS_900_a, 0.0, endTem_900), (YS_900_b, 0.1, endTem_900), ))
mdb.models[modelname].materials['Ti'].Expansion(temperatureDependency=ON, 
    table=((alphai_300, endTem_300), (alphai_325, endTem_325), (alphai_350, endTem_350), (alphai_375, endTem_375), 
    (alphai_400, endTem_400), (alphai_425, endTem_425), (alphai_450, endTem_450), (alphai_475, endTem_475), 
    (alphai_500, endTem_500), (alphai_525, endTem_525),  (alphai_550, endTem_550),  (alphai_575, endTem_575), 
    (alphai_600, endTem_600), (alphai_625, endTem_625), (alphai_650, endTem_650), (alphai_675, endTem_675), 
    (alphai_700, endTem_700), (alphai_725, endTem_725), (alphai_750, endTem_750), (alphai_775, endTem_775), 
    (alphai_800, endTem_800), (alphai_825, endTem_825),  (alphai_850, endTem_850),  (alphai_875, endTem_875), 
    (alphai_900, endTem_900)), zero=298.15)
mdb.models[modelname].materials['Ti'].Conductivity(temperatureDependency=ON,
     table=((Ki_300, endTem_300), (Ki_325, endTem_325), (Ki_350, endTem_350), (Ki_375, endTem_375), 
     (Ki_400, endTem_400), (Ki_425, endTem_425), (Ki_450, endTem_450), (Ki_475, endTem_475), 
     (Ki_500, endTem_500), (Ki_525, endTem_525), (Ki_550, endTem_550), (Ki_575, endTem_575), 
     (Ki_600, endTem_600), (Ki_625, endTem_625), (Ki_650, endTem_650), (Ki_675, endTem_675), 
     (Ki_700, endTem_700), (Ki_725, endTem_725), (Ki_750, endTem_750), (Ki_775, endTem_775), 
     (Ki_800, endTem_800), (Ki_825, endTem_825), (Ki_850, endTem_850), (Ki_875, endTem_875), 
     (Ki_900, endTem_900)))
     
     
######�����������######
#sectionLayer1��sectionLayer(cengshu)�����ڲ㵽���
d = OrderedDict()
for i in range(1, cengshu+1):
    matname='O-I-'+str(i)
    d["sectionLayer{0}".format(i)] = section.SectionLayer(material=matname, thickness=CH, orientAngle=0.0, numIntPts=InPt, plyName='')


laylist=[]
for key, val in d.items():
    laylist.append(d[key])


m.CompositeShellSection(name='Shell_shell_wgjj-'+ptname, preIntegrate=OFF, 
    idealization=NO_IDEALIZATION, symmetric=False, thicknessType=UNIFORM, 
    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=POINTWISE, 
    nTemp=InPt, useDensity=OFF, integrationRule=SIMPSON, layup=(laylist))
            
m.HomogeneousShellSection(thickness=tw,material='Ti',name='Shell_jintiao_wgjj-'+ptname,
    nodalThicknessField='',numIntPts=5,integrationRule=SIMPSON,
    poissonDefinition=DEFAULT,preIntegrate=OFF,temperature=GRADIENT,
    thicknessModulus=None,useDensity=OFF)                           #��������   #ptname='p1'

######������ͼ######    
s=m.ConstrainedSketch(name='rib1',sheetSize=4*pi*R)
s1=m.ConstrainedSketch(name='rib2',sheetSize=4*pi*R)
s2=m.ConstrainedSketch(name='shell',sheetSize=4*pi*R)

######��������(����)������ʵ��######
pt1=(R,0.0)
pt2=(R-h,0.0)
    
s.Line(point1=pt1,point2=pt2)
s.Spot(point=(0.0,0.0))
s.radialPattern(geomList=(s.geometry[2], ),vertexList=(),number=vn,
    totalAngle=360.0,centerPoint=(0.0,0.0))  #��s��ͼ�н����ߣ��ڻ������У�����Ϊvn
        
p1=m.Part(name='rib1_temp_temp',dimensionality=THREE_D,type=DEFORMABLE_BODY)
p1.BaseShellExtrude(sketch=s,depth=H)   #�ɲ�ͼs��������p1���߶�H
facelist99=[]
for face in p1.faces:
    facelist99.append(p1.faces.findAt(face.pointOn))
p1.Set(name='Set_Face_for_section_ce1',faces=facelist99)  #��p1��������ļ���
r.Instance(name='rib1_merge_merge',part=p1,dependent=ON)  #��p1���ɷǶ���ʵ��
    
s1.CircleByCenterPerimeter(point1=pt1,center=(0.0,0.0))
s1.CircleByCenterPerimeter(point1=pt2,center=(0.0,0.0))  #��s1��ͼ������һ��Բ��
p3=m.Part(name='rib3_temp_temp',dimensionality=THREE_D,type=DEFORMABLE_BODY)
p3.BaseShell(sketch=s1)          #���ڲ�ͼs1����һ��3άʵ��ǵ�Ԫ����p3
facelist99=[]
for face in p3.faces:
    facelist99.append(p3.faces.findAt(face.pointOn))
p3.Set(name='Set_Face_for_section_ce2',faces=facelist99)   #��p3��������һ������Set_Face_for_section_ce2
r.Instance(name='rib3_merge_merge',part=p3,dependent=ON)   #��p3����һ���Ƕ���ʵ��
r.LinearInstancePattern(instanceList=('rib3_merge_merge', ),
    direction1=(0.0,0.0,10.0),
    direction2=(0.0,1.0,0.0),number1=hn,number2=1,
    spacing1=hgap,spacing2=2000.0)       #��ʵ��p3��z�������У���ĿΪhN������������Ŀ��ǰ���Ѷ��壩��ƫ�ƾ���hgap������2Ϊy�򣬲�����
    
    
#���������ɵĲ���ʵ�кϲ����ɲ���ʵ�У�ԭ�����ɵ�ʵ��ɾ������������ɾ��������Ȼͬʱ������Ӧ������
r.InstanceFromBooleanMerge(name='Jintiao',instances=
    r.instances.values(),originalInstances=DELETE)
p=m.parts['Jintiao']            #(����)��������
p.Set(name='Set_Face_for_section_rib_wgjj',faces=p.faces)
p.SectionAssignment(region=p.sets['Set_Face_for_section_rib_wgjj'],sectionName='Shell_jintiao_wgjj-'+ptname) 
jt='Jintiao-1'




######��������(��Ƥ)������ʵ��###### 
m.ConstrainedSketch(name='__profile__',sheetSize=2000.0)
m.sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0,0.0),point1=(R,0.0))
m.Part(dimensionality=THREE_D,name='SHELL-G',type=
    DEFORMABLE_BODY)
m.parts['SHELL-G'].BaseShellExtrude(depth=H,sketch=
    m.sketches['__profile__'])
del m.sketches['__profile__']
p=m.parts['SHELL-G']         #��Ƥ����
p.Set(name='Set_Face_for_section_shell_wgjj',faces=p.faces)
       
m.rootAssembly.DatumCsysByDefault(CARTESIAN)    #����Ĭ������ϵΪ�翨������ϵ
m.rootAssembly.Instance(dependent=ON,name='SHELL-G-1',part=m.parts['SHELL-G'])  #��Ƥʵ��



gg='SHELL-G-1'     #��Ƥʵ��
ff='SHELL-G'

r=m.rootAssembly                                                                                
r.InstanceFromBooleanMerge(name='Shell-allall',instances=(r.instances[ff+'-1'],
    r.instances[jt], ),keepIntersections=ON,originalInstances=DELETE,domain=GEOMETRY)
     
ff='Shell-allall'
gg=ff+str(-1)
p=m.parts[ff]


######������R�����棨�ϲ�����Ƥ�ϵ��棩######
tol=1e-2              #������Χ
facelist=[]           #��������б�

for face in p.faces:
    tmp=True                                           #�жϱ���tmp
    for j in face.getVertices():                       #�����и����������j��
        pt=p.vertices[j].pointOn                       #��vertices�б���ĵ�Ը���pt
        d2=pow(pt[0][0],2)+pow(pt[0][1],2)             #pow��x��y��z������x��y�η��ڳ���z���������pow(pt[0][0],2)����pt��x��ƽ����
        if abs(d2-R*R)>tol:                            #abs()��ȡ����ֵ�����d2��pt��ԭ�㳤�ȵ�ƽ����ȥr��ƽ����ֵ����tol������������r�ķ�Χ�⡣
            tmp=False                                  #tmpΪ����
    if tmp==True:                                      #���������R��Χ�࣬�⽫���ϵõ��ĵ�����facelist�С�
        facelist.append(p.faces.findAt(face.pointOn))  #��������FindAt����ͨ��faces �ϵ�����ķ�����λ����Ӧ���棬���������б��С�findAtΪ��λ����ͨ��ĳһ������λ�ҵ��������
p.Set(name='Set_Face_shell_wgjj',faces=facelist)       #��facelist�б���Ϊ����Set_Face_shell_wgjj�� p=m.parts[ff]

p.SectionAssignment(region=p.sets['Set_Face_shell_wgjj'], sectionName='Shell_shell_wgjj-'+ptname, offset=0.0, 
    offsetType=BOTTOM_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


######�������������ᣬɸѡ���ǵ׶˺Ͷ��˱�Ե######
c_name='z'
c_value1=0.0       #������Ӧ��������ϵ��������Χ��������Χ��c_value1��c_value2����
c_value2=None
CordType='REC'        #ֱ������ϵ

tol=1.0e-5

if c_value2==None:
    c_value2=c_value1

if c_name=='x':
    index=0
elif c_name=='y':
    index=1
elif c_name=='z':
    index=2

edgelist=[]

for i in range(0,p.edges[-1].index+1):
    tmp=True
    for j in p.edges[i].getVertices():
        pt=p.vertices[j].pointOn
        if CordType=='REC':
            value=pt[0][index]
        
        if CordType=='CYL':
            if index==0:
                value=sqrt(pt[0][0]**2+pt[0][1]**2)
            if index==1:
                value=atan(pt[0][1]/pt[0][0])
                
                
        if value<(min(c_value1,c_value2)-tol) or value>(max(c_value1,c_value2)+tol):
                tmp=False
    if tmp:
        edgelist.append(p.edges.findAt(p.edges[i].pointOn))      

p.Set(name='Set_Edge_shellbot_wgjj',edges=edgelist)   #�׶˵㼯��   

c_name='z'
c_value1=H
c_value2=None
CordType='REC'

tol=1.0e-5

if c_value2==None:
    c_value2=c_value1

if c_name=='x':
    index=0
elif c_name=='y':
    index=1
elif c_name=='z':
    index=2

    
edgelist=[]

for i in range(0,p.edges[-1].index+1):
    tmp=True
    for j in p.edges[i].getVertices():
        pt=p.vertices[j].pointOn
        if CordType=='REC':
            value=pt[0][index]
        
        if CordType=='CYL':
            if index==0:
                value=sqrt(pt[0][0]**2+pt[0][1]**2)
            if index==1:
                value=atan(pt[0][1]/pt[0][0])
                
                
        if value<(min(c_value1,c_value2)-tol) or value>(
            max(c_value1,c_value2)+tol):
                tmp=False
    if tmp:
        edgelist.append(p.edges.findAt(p.edges[i].pointOn))     

p.Set(name='Set_Edge_shelltop_wgjj',edges=edgelist)    #���˵㼯��


######�������������ᣬɸѡ���ǽ���######
#ɸѡ����R��R-h������������Χ��ı�Ե �������嵽����Set_Edge_midedge_wgjj       
p=m.parts[ff]    #����Ƥ��jintiao���ɵĲ�������p
tol=1.0e-5       #�趨������ѡ��
edgelist=[]      #�����б�edgelist

for i in range(0,p.edges[-1].index+1):      #��p�и�����Ե�������ŵ���ABAQUS�У�����˳ʱ���ţ�������ʱ�룬�����ɵ��ȱ�ţ�
    tmp=True                                #����tmpĬ��Ϊture
    for j in p.edges[i].getVertices():       #���õ��ı�Ե��������ѭ������j��
        pt=p.vertices[j].pointOn             #����x2+y2�Ŀ�������Ϊvalue�����ж�
        value=sqrt(pt[0][0]**2+pt[0][1]**2)
                     
        if value>R-h-tol and value<R+tol:          #����and�жϷ��Ž����ж�
            tmp=False
        if tmp==False:                    #���ڱ�Ե�ϵĵ㵽ԭ�����ΪR��R-h��ĵ�����б��У�����ֹ���ҵ��ѭ��������һ����ѭ������Ѱ�ұ�Ե�����ս��ҵ��ı�Ե�����б���һ������Set_Edge_midedge_wgjj
            edgelist.append(p.edges.findAt(p.edges[i].pointOn))
            break                  #break����ֹѭ����
p.Set(name='Set_Edge_midedge_wgjj',edges=edgelist)

    
#��������Ե�ϵĵ��������ڰ뾶R-h�ϵ��ҵ���������Ϊ�ڱ�Ե���ϣ�   
c_name='x'       #Ĭ����������ϵ����x��
c_value1=R-h      #���ֵΪR-h
c_value2=None
CordType='CYL'      #����������ϵ

tol=1.0e-5         #�������

if c_value2==None:      #���û�и���c_value2ֵ����c_value2=c_value1=R-h
    c_value2=c_value1
if c_name=='x':      #��x��y��z����ϵ��0,1,2(����index)����
    index=0
if c_name=='y':
    index=1
if c_name=='z':
    index=2
    
edgelist=[]             #�����Ե�б�

for i in range(0,p.edges[-1].index+1):    #����Ե����ѭ������i��
    tmp=True                             
    for j in p.edges[i].getVertices():  #����Ե�Ķ�������ѭ������j������j�����������pt
        pt=p.vertices[j].pointOn        #����ֱ������ϵ�������궨�岻ͬ��valueֵ�Ա����ɸѡ
        if CordType=='REC':
            value=pt[0][index]
        
        if CordType=='CYL':               #��������õ��Ǿ��루c_nameΪx��������ǣ�c_nameΪy��
            if index==0:
                value=sqrt(pt[0][0]**2+pt[0][1]**2)
            if index==1:
                value=atan(pt[0][1]/pt[0][0])
                
                
        if value<(min(c_value1,c_value2)-tol) or value>(
            max(c_value1,c_value2)+tol):
                tmp=False                   #������c_value1��c_value2֮����tmpΪFalse���������
    if tmp:                               #�����c_value1��c_value2��Χ�࣬����Ӧ��Ե��������б��У������弯�� Set_Edge_inedge_wgjj���ڱ�Ե��
        edgelist.append(p.edges.findAt(p.edges[i].pointOn))
p.Set(name='Set_Edge_inedge_wgjj',edges=edgelist)

#ɸѡ���ڰ뾶R�ϵı�Ե��������Ϊ���Ե��
c_name='x'         #ͬ��Ĭ����������ϵΪx�ᣬ������ʼֵΪc_value1=R��
c_value1=R
c_value2=None
CordType='CYL'    #Ĭ������ϵΪ������ϵ

tol=1.0e-5        #�������

if c_value2==None:     ##���û�и���c_value2ֵ����c_value2=c_value1=R-h���Լ���x��y��z������ֵ0,1,2����
    c_value2=c_value1
if c_name=='x':
    index=0
if c_name=='y':
    index=1
if c_name=='z':
    index=2
    
edgelist=[]                 #����һ���ձ�Ե�б�

for i in range(0,p.edges[-1].index+1):  #ͬ������ѭ��p����Ƥ+JINGTIAO�������еı�Ե�������ݱ�Ե�ϵĸ�������Ĳ�ͬ����ϵ�µ�Valueֵ����ɸѡ
    tmp=True
    for j in p.edges[i].getVertices():
        pt=p.vertices[j].pointOn
        if CordType=='REC':
            value=pt[0][index]
        
        if CordType=='CYL':
            if index==0:
                value=sqrt(pt[0][0]**2+pt[0][1]**2)
            if index==1:
                value=atan(pt[0][1]/pt[0][0])
                
                
        if value<(min(c_value1,c_value2)-tol) or value>(
            max(c_value1,c_value2)+tol):
                tmp=False                 #����Щ��Ե����R����Χ��ı�Ե�ų�
    if tmp:                               #��Ե�ϵ���R����Χ�ڵı�Ե�����б������ս��б��е�ֵ���뼯��Set_Edge_outedge_wgjj�����Ե��
        edgelist.append(p.edges.findAt(p.edges[i].pointOn))
p.Set(name='Set_Edge_outedge_wgjj',edges=edgelist)


######���֣���������######
#���񻮷֣��ȶ��岿�����ӣ����ڽ��н������Ĳ��������������Ե���Ͻ����� ��ͬ��Ϊ�߲��֣����б�Ե���Ͻ��������Ϊ�߲��֣���������Լ�����á��������񻮷��ҽ��б��档
p.seedPart(size=msize1,deviationFactor=0.1)  #���廮�������ȫ�����ӳߴ�Ϊmsize1�����ƫ������Ϊ0.1
p.setMeshControls(regions=p.faces,elemShape=QUAD_DOMINATED,technique=FREE,allowMapped=True) #�������Կ��ƣ���Ԫ��״�ı���Ϊ�������ּ������ɻ��֣������㷨�Լ��ں��ʵط�����ʹ��ӳ������

p.seedEdgeBySize(edges=p.sets['Set_Edge_midedge_wgjj'].edges,size=msize2,deviationFactor=0.1,constraint=FIXED)  #Ϊ�ߣ�����Ե�����֣�3�������ƫ�����Ӷ�Ϊ0.1
p.seedEdgeBySize(edges=p.sets['Set_Edge_outedge_wgjj'].edges,size=msize1,deviationFactor=0.1,constraint=FIXED)  #���Ե���ڱ�Ե���Ƶ�Ԫ�ߴ�Ϊmsize1�������Ե�ߴ�Ϊmsize2
p.seedEdgeBySize(edges=p.sets['Set_Edge_inedge_wgjj'].edges,size=msize1,deviationFactor=0.1,constraint=FIXED)   #Լ����������Լ������Fixed����ִ�ģ��̶�����ģ���������ı䵥Ԫ��
    
p.generateMesh(seedConstraintOverride=ON)  #�������������л��֣������Ӳ�����Լ�����¸��ǣ������أ���

#mdb.saveAs(pathName='aaa')    #����ģ���ļ���·�����ƣ�aaa

#######���òο���#######
r.ReferencePoint(point=(0.0,0.0,H))
r1 = r.referencePoints
refPoints1=r1.keys()[-1]
key1=r1[refPoints1]
r.Set(referencePoints=(key1,), name='Set_Rpt_wgjj')   #�˲��������

m.Coupling(name='Coupling_top_wgjj', controlPoint=r.sets['Set_Rpt_wgjj'], 
    surface=r.instances[gg].sets['Set_Edge_shelltop_wgjj'], 
    influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)  #�˲���� 

############################################
#mdb.saveAs(pathName='linear_buckling')
    
nodelist=p.nodes
p.Set(name='Set_Face_nodelist_wgjj',nodes=nodelist)

facelist=p.faces
edgelist=p.edges
verticelist=p.vertices
p.Set(name='Set_shell_allall',faces=facelist,edges=edgelist,vertices=verticelist)
  
# mdb.saveAs(pathName='linear_buckling') 


######����������######
if computmethod=='Explicit':        #��ʾ����ѧ���� 
    loadValue=-H*0.01
    m.ExplicitDynamicsStep(name='Step-1',previous='Initial',timePeriod=computTime)    #����������Step-1
    m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=40,
        timeMarks=ON)
    m.TabularAmplitude(name='Amp-1',timeSpan=STEP,
        smooth=SOLVER_DEFAULT,data=((0.0,0.0),(computTime,1.0)))   #�ڹ���������һ������Ϊ��ķ�ֵ����
    m.DisplacementBC(createStepName='Step-1',fieldName='',
        fixed=OFF,localCsys=None,name='disp',region=r.sets['Set_Rpt_wgjj'],u1=UNSET,
        u2=UNSET,u3=loadValue,ur1=UNSET,ur2=UNSET,ur3=UNSET,amplitude='Amp-1')      #��Step-1��ʩ�Ӻ���(������ϵ�)
elif computmethod=='Linear_Buck':
    m.StaticStep(name='Step-1', previous='Initial')  #, initialInc=0.05, nlgeom=ON
    m.StaticStep(name='Step-2', previous='Step-1')
    m.BuckleStep(name='Step-3', previous='Step-2', numEigen=1, 
        eigensolver=LANCZOS, minEigen=1e4, blockSize=DEFAULT, maxBlocks=DEFAULT)   #�趨һ����С�½� ���������  None  
    m.Temperature(name='Predefined Field-1', 
        createStepName='Initial', region=r.instances[gg].sets['Set_shell_allall'], distributionType=UNIFORM, 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Tem, ))     #��ʼ�¶ȳ�
    m.Temperature(name='Predefined Field-2', 
        createStepName='Step-1', region=r.instances[gg].sets['Set_Face_for_section_rib_wgjj'], 
        distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Temshell[TemPt-1], ))  #���������¶ȳ�
    m.Temperature(name='Predefined Field-3', 
        createStepName='Step-1', region=r.instances[gg].sets['Set_Face_for_section_shell_wgjj'], 
        distributionType=UNIFORM, crossSectionDistribution=POINTS_THROUGH_SECTION, magnitudes=(Temlist))   #��Ƥ�����¶ȳ�
    m.ConcentratedForce(name='Load-1', createStepName='Step-2', 
        region=r.sets['Set_Rpt_wgjj'], cf3=-Jzhz, distributionType=UNIFORM, field='', localCsys=None)      #��Ч�¶Ⱥ���       
    m.ConcentratedForce(name='Load-2', createStepName='Step-3', 
        region=r.sets['Set_Rpt_wgjj'], cf3=-1.0, distributionType=UNIFORM, field='', localCsys=None)       #������������
elif computmethod=='Nonlinear_Newton':
    m.StaticStep(name='Step-1', previous='Initial',nlgeom=ON)
    m.Temperature(name='Predefined Field-1', 
        createStepName='Initial', region=r.instances[gg].sets['Set_shell_allall'], distributionType=UNIFORM, 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Tem, ))    #��ʼ�¶ȳ�
    m.Temperature(name='Predefined Field-2', 
        createStepName='Step-1', region=r.instances[gg].sets['Set_Face_for_section_rib_wgjj'], 
        distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Temshell[TemPt-1], ))  #���������¶ȳ�
    m.Temperature(name='Predefined Field-3', 
        createStepName='Step-1', region=r.instances[gg].sets['Set_Face_for_section_shell_wgjj'], 
        distributionType=UNIFORM, crossSectionDistribution=POINTS_THROUGH_SECTION, magnitudes=(Temlist))    #��Ƥ�����¶ȳ�

    
m.DisplacementBC(amplitude=UNSET,createStepName='Initial',distributionType=UNIFORM,
    fieldName='',fixed=OFF,localCsys=None,name='BC-1',region=r.instances[gg].sets['Set_Edge_shellbot_wgjj'],
    u1=0.0,u2=0.0,u3=0.0,ur1=0.0,ur2=0.0,ur3=0.0)              #���Լ��(�ײ�)   #Լ������ص�����ϵҪһ��
m.DisplacementBC(amplitude=UNSET,createStepName='Initial',distributionType=UNIFORM,
    fieldName='',fixed=OFF,localCsys=None,name='BC-2',region=r.sets['Set_Rpt_wgjj'],
    u1=0.0,u2=0.0,u3=UNSET,ur1=0.0,ur2=0.0,ur3=0.0)            #���Լ��(����)

  
######ɾ������Ĳ�����ģ��######    
for key in m.parts.keys():
    if key.find('temp_temp')!=-1:   #��temp_temp���ǲ������ֵ�ĩβ����ɾ���ò���
        del m.parts[key]


for key in mdb.models.keys():
    if mdb.models[key].parts.keys()==[]:
        del mdb.models[key]

       
m.rootAssembly.regenerate()    #ɾ�������ģ�͸�װ������������һ�Σ���������Ч�� 
m.rootAssembly.regenerate()


######ģ���������######
zd=mdb.models[modelname].rootAssembly.getMassProperties()
weight=zd['mass']
weight2=weight*1000.0
outfile=open('output_W.txt','w')
outfile.write('wt='+str(weight2)+'\n')
outfile.close()

#######����Job����######
mdb.Job(name='am22', model=modelname, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=10, 
    numDomains=20, numGPUs=0)
mdb.saveAs(pathName='linear_buckling2')
    
mdb.jobs['am22'].submit(consistencyChecking=OFF)

mdb.jobs['am22'].waitForCompletion()


# �Ե�ǰ�ӿ��е�������ݿ���в���
odb=session.openOdb(name='E:/HSM_AK_FGMshell/SSM_lightweight/am22.odb')

# �������ģ̬����������ֵ
lastFrame=odb.steps['Step-3'].frames[1].description     
            
outfile=open('output2_F.txt','w')
outfile.write('lingjieli='+str(lastFrame)+'\n')
outfile.close()
