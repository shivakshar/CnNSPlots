

import numpy as np
import math
#from matplotlib import interactive
#interactive(True)
import matplotlib.pyplot as plt
import pandas as pd
import sys as syst
from CNnsCS import *

g_b=pd.read_csv('./gamma_bkg.txt', sep='\t', names=["bin_center","# events"],index_col="bin_center")
n_b=pd.read_csv('./neutron_bkg.txt', sep='\t', names=["bin_center","# events"],index_col="bin_center")

e_w=(g_b.index[1]-g_b.index[0])*1e2*5
e_l=g_b.index*1e2*5-(e_w/2)
e_r=g_b.index*1e2*5+(e_w/2)
e_l=list(e_l)
e_l.append(e_r[-1])
g_n_b=(g_b.values+n_b.values).transpose()[0]

g=g_b.values.transpose()[0]
n=n_b.values.transpose()[0]

def num_b(t1,t2):#eV 
    for i in range(len(e_l)):
        if e_l[i]>t1:
            a=i-1
            break
        

    for i in range(len(e_l)):
        if e_l[i]>t2:
            b=i
            break
    
    
    #print(a,b)    
    #print(e_l[a-1],e_l[b-1])
    
    lo=g_n_b[a]*(e_l[a+1]-t1)/e_w
    mid=sum(g_n_b[a+1:b])
    hi=g_n_b[b-1]*(t2-e_l[b])/e_w
    #print(lo)
    #print(mid)
    #print(hi)
    
    return lo+mid+hi
    

def chicha(bins,s_b,sys,t,th):#th in eVnr; Rest of code in keV
#    s_total=rToth[th]
#    b_total=(s_total/s_b)*t
#    b=[b_total/bins]*bins
    #print("Eval Integral")    
    s=[C_T1_T2(i,i+((1-(th/1000))/bins))*(60*60*24)*t*m/(L*L) for i in np.linspace((th/1000), 1, num=bins)]
    #print(s)
    #print(b)
    b=np.array([num_b(i,i+((1000-(th))/bins))*(60*60*24)*t*m/(L*L) for i in np.linspace((th), 1000, num=bins)])
    b=(b/sum(b))*(sum(s)/s_b)
    
    b_sig=np.vectorize(lambda b,sys: (b+((sys*b)**2))**0.5)
    b_sig_val=b_sig(b,sys)
    #print(b_sig_val)
    chisq_val=np.vectorize(lambda sys:(s/b_sig_val)**2)
    return(sum(chisq_val(sys)))

chicha=np.vectorize(chicha)

def chi_plot(th,s_b,sys,bins):
    t=10**np.linspace(-2,3,10)
    chis=chicha(bins,s_b,sys,t,th)
    plt.plot(t,chis)
    plt.xscale("log")
    plt.xlabel("Time (days)", size=20)
    plt.ylabel(r"$\Delta\chi^2$", size=20)
    plt.ylim(0,30)
    plt.plot((t[0],t[-1]),(9,9),linestyle="--",c='k')
    plt.plot((t[0],t[-1]),(25,25),linestyle="--",c='k')
    plt.title("Comparing Signal and Background distributions",size=15)
    plt.tick_params(labelsize=15)
    plt.text(200,9,r"$3\sigma$",size=25)
    plt.text(200,25,r"$5\sigma$",size=25)
    plt.title(details)
    plt.tight_layout()
    plt.savefig('./Plots/chi2plot.png')
    #plt.show(block=False)
    #plt.ion()

th=float(syst.argv[1])
s_b=float(syst.argv[2])
sys=float(syst.argv[3])
bins=float(syst.argv[4])
L=float(syst.argv[5])
m=float(syst.argv[6])
details="distance to core="+str(L)+"m| mass of det="+str(m)+"kg\nthreshold="+str(th)+"eVnr| S/B="+str(s_b)+"\nsyst err="+str(sys*1e2)+"%| # of bins="+str(bins)

chi_plot(th,s_b,sys,bins)



