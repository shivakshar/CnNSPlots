#Probability at L=1m
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chisquare
import sys as syst
#from CnNS import pdf_E
#from CnNS import random_1
#from CnNS import pdf_slice_1  

from CNnsCS import *
#from CnNS_Counting_2 import *


#plt.plot(pdf_E,random_1)
#plt.plot(pdf_E,pdf_slice_1) 

# E in eV->keV
#L=1#m
#bins=100
#th=10#eV
t=300#day

th=float(syst.argv[1])
s_b=float(syst.argv[2])
sys=float(syst.argv[3])
bins=float(syst.argv[4])
L=float(syst.argv[5])
m=float(syst.argv[6])
details="distance to core="+str(L)+"m| mass of det="+str(m)+"kg\nthreshold="+str(th)+"eVnr| S/B="+str(s_b)+"\nsyst err="+str(sys*1e2)+"%| # of bins="+str(bins)
#print(len(pdf_E))
#print(len(random_1))
#print(len(pdf_slice_1))

s=[C_T1_T2(i,i+((1-(th/1000))/bins))*(60*60*24)*t*m/(L*L) for i in np.linspace((th/1000), 1, num=bins+1)[:-1]]
pdf_E=np.linspace((th/1000), 1, num=bins+1)[:-1]
pdf_slice_1=random_1=s

def Oscillate(pdf_E,pdf_slice_1,a,m2):
    k=1.27 #Constant associated with Units   
    P=a*(np.sin(m2*k*((L/1e3)/(pdf_E/1e6)))**2)#length from m to km, energy from keV to GeV
    pdf_1_osc=pdf_slice_1*(1-P)#Energy Flux with Oscillations
    return pdf_1_osc

#plt.plot(pdf_E,pdf_slice_1); plt.plot(pdf_E,Oscillate(pdf_E,pdf_slice_1,0.3,1e-2))


adelta=0.01
m2delta=0.01
a=np.arange(0+adelta, 1+adelta, adelta)
m2=np.arange(0+m2delta, 10+m2delta, m2delta)

A,M2=np.meshgrid(a, m2)

sys=0
def calc_chi_sq(Obs,Exp,sys):
    Exp_sig=np.vectorize(lambda Exp,sys: (Exp+((sys*Exp)**2))**0.5)
    Exp_sig_val=Exp_sig(Exp,sys)
    #print(b_sig_val)
    chisq_val=np.vectorize(lambda sys:((Obs-Exp)/Exp_sig_val)**2)
    return(sum(chisq_val(sys)))
    

def Chi(pdf_E,pdf_slice_1,a,m2,random_1):
    v=calc_chi_sq(Oscillate(pdf_E,pdf_slice_1,a,m2),random_1,sys)
    #v=chisquare(Oscillate(pdf_E,pdf_slice_1,a,m2),f_exp=random_1,ddof=len(pdf_E))
    return (v)

Data=[]
for ax in a:
    temp=[]
    for m2x in m2:
        temp.append(Chi(pdf_E,pdf_slice_1,ax,m2x,random_1))
    Data.append(temp)

CHI=np.transpose(np.array(Data))

lines=[]
labels=[]

plt.figure(figsize=(16,16))


cs=plt.contour(A,M2,CHI, levels=[9,25], colors=['r','b']) #Specific confidence levels correspond to specific chi-sq values
lines=[cs.collections[0],cs.collections[1]]
labels = ['3 sigma','5 sigma']
plt.legend(lines, labels, loc=3)

#plt.pcolor(A,M2,CHI)
#plt.colorbar()
plt.xlabel(r"sin$^2(2 \theta)$",size=20)
plt.ylabel(r"$\Delta m^2 \ ($eV$^2)$",size=20)
plt.xlim(0.01,1)
plt.ylim(0.01,10)
plt.xscale("log")
plt.yscale("log")
plt.tick_params(labelsize=20)
plt.title('Runtime of '+str(t)+' days\n'+details,size=10)#if needed, change Title Appropriately
#plt.figure(figsize=(16,12))
plt.tight_layout()
plt.savefig("./Plots/SterileSensitivity.png")