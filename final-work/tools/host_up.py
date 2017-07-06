#!/usr/bin/python
import sys
from subprocess import call
#labos={6:'10.2.6.1-25',7:'10.2.6.1-25'}
#f=open('.nmap_out.txt','w')
#call(['nmap','-sP labos[int(sys.argv[1])]'],stdout=f)
#f.close()
#f=open('.nmap_out.txt','r')
f=open(sys.argv[1],'r')
i=0
lineas=""
for l in f.readlines():
    if i%2==0:
        lineas=lineas+l
    i=i+1
f.close()
#call(['rm','.nmap_out.txt']) 
lineas2=lineas.split("(")
lineas3=[]
for i in range(len(lineas2)):
        lineas3.append((lineas2[i]).split(")")[0])
ip=lineas3[1:len(lineas3)-1]
f=open(sys.argv[2],'w')
for elt in ip:
    f.write(elt+"\tslots=2\n")
f.close()    
