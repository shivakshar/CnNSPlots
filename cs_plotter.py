import numpy as np
import math
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

def CSpres_2(S_B,sys,bins,th,n,t):#t in days
    #s_total=rToth[th]
    #print(s_total)
    #b_total=(s_total/S_B)*t
    #print(b_total)
    
    s_b=S_B
    #print("Eval Integral")    
    sig=np.vectorize(lambda x,sy: (x+(x*x*sy*sy))**0.5) 
    s=np.array([C_T1_T2(i,i+((1-(th/1000))/bins))*(60*60*24)*t for i in np.linspace((th/1000), 1, num=bins)])
    
    b=np.array([num_b(i,i+((1000-(th))/bins))*(60*60*24)*t for i in np.linspace((th), 1000, num=bins)])
    b=(b/sum(b))*(sum(s)/s_b)
    #print(s)
    #print(b)
    #s_to_b=s/b
    #k=b/s
    #print(s_to_b)
    
    
    #print(sig(s,sys))
    #print(sig(b,sys))
    #mu_obs_2=(((sig(s,sys)*sig(s,sys)*s)+(sig(b,sys)*sig(b,sys)*(((s+b)-(n*sig(s+b,sys)))-b)))/((sig(s,sys)*sig(s,sys))+(sig(b,sys)*sig(b,sys))))**2
    sig_obs_2=((sig(s,sys)*sig(s,sys)*sig(b,sys)*sig(b,sys))/((sig(s,sys)*sig(s,sys))+(sig(b,sys)*sig(b,sys))))
    #print(sig_obs_2)    
    mu_obs_2=(((sig(s,sys)*sig(s,sys)*s)+(sig(b,sys)*sig(b,sys)*(((s+b)+(n*(sig(s+b,sys))))-b)))/((sig(s,sys)*sig(s,sys))+(sig(b,sys)*sig(b,sys))))**2
    #print(s)
    #print(b)    
    #print(s+b)
    #print(sig(s+b,sys))
    #print((s+b)+(n*(sig(s+b,sys))))
    #print((((s+b)+(n*(sig(s+b,sys))))-b))
    
    prec=(sig_obs_2/mu_obs_2)
    
    #print(1/prec)
    #print((((1/sum(1/prec))+(sys*sys))**0.5))    
    return(((1/sum(1/prec))+(sys*sys))**0.5)#Add the systematic error


def CSpres_plot_2(S_B,sys,bins,th):#eV
    t=10**np.linspace(-1,3,10)
    #min_CS_pres=CSpres_2(S_B,sys,bins,th,-5,t)
    mean_CS_pres=CSpres_2(S_B,sys,bins,th,0,t)
    #max_CS_pres=CSpres_2(S_B,sys,bins,th,+5,t)
    
    #plt.fill_between(t,min_CS_pres,max_CS_pres, alpha=0.1, color="b")
    plt.plot(t,mean_CS_pres,linestyle="--", color='b')
    plt.plot(t,[sys]*len(t),c="r")
    plt.ylim(0,1)
    #plt.yscale("log")
    plt.xscale("log")    
    plt.xlabel("Time (days)")
    plt.ylabel(r"$(\frac{d\sigma}{\sigma})$", size=20)
    #print(min_CS_pres)
    #print(mean_CS_pres)
    #print(max_CS_pres)
    plt.title(details)
    plt.tight_layout()
    plt.savefig('./Plots/crosssectionplot.png')
    #plt.show()
    #syst.exit()

th=float(syst.argv[1])
s_b=float(syst.argv[2])
sys=float(syst.argv[3])
bins=float(syst.argv[4])
L=float(syst.argv[5])
m=float(syst.argv[6])

details="distance to core="+str(L)+"m| mass of det="+str(m)+"kg\nthreshold="+str(th)+"eVnr| S/B="+str(s_b)+"\nsyst err="+str(sys*1e2)+"%| # of bins="+str(bins)
#print("H")
CSpres_plot_2(s_b,sys,bins,th)

