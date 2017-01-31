import numpy as np
import math
import matplotlib.pyplot as plt


##Input Parameters for Sterile Neutrino Sensitivity Study

L=1              #meters from reactor core
t=1          #number of years
#N_data=100     #Number of data points generated()
th=0.01           #detector threshold in keV

#------------------------------
"""
tyr=60*60*24*365
t=t*tyr
"""
s2z=0.23120
pNC_vN=1.0086
k_vN=0.9978
Lul=-0.0031
Ldl=-0.0025
LdR=7.5e-5
LuR=0.5*(7.5e-5)

Z=32
N=40

gp_V=pNC_vN*(0.5-(2*k_vN*s2z))+(2*Lul)+(2*LuR)+Ldl+LdR
gn_V=-(0.5*pNC_vN)+Lul+LuR+(2*Ldl)+(2*LdR)


Nplus_Nminus=2
Zplus_Zminus=-2


G_V=(gp_V*Z)+(gn_V*N)
G_A=(gp_V*(Zplus_Zminus))+(gn_V*(Nplus_Nminus))

M=66.9792e6 #keV: Mass of Nucleus
M_det=5.60959e32 #keV: Mass of Detector 1kg


#NSI parameters: 0

from scipy.integrate import quad
from scipy.stats import maxwell

flux_1m=1.511067123220e12#/cm^2/s at 1m
flux_r=lambda r: flux_1m/(r*r) #r in meters
total_flux=flux_r(L)#/cm^2/s at 1m

#np.arange(1e2,20e3,1e2)#E=100keV - 20MeV with 100keV step 

form_n_flux=lambda En:maxwell.pdf((En+1.9e3)/1.7e3)#En in keV
total=quad(lambda En:form_n_flux(En),0,np.inf)[0] 

norm_n_flux=lambda En:form_n_flux(En)/total
total_norm_n_flux=quad(norm_n_flux,0,np.inf)[0]

n_flux=lambda En: total_flux*norm_n_flux(En)
total_n=quad(n_flux,0,np.inf)[0]



#t=1#second
Gf2=(1.16637e-17)**2#keV-4
#area_det=34.322512#cm^2  /((4)*math.pi*100*100)

C_T=lambda En,T: (Gf2*M/(2*math.pi)) * ( ((G_V+G_A)**2) + (((G_V-G_A)**2)*((1-(T/En))**2)) - (((G_V**2)-(G_A**2))*(M*T/(En*En))) )

nat_unit_fac=(1.97e-8)**2#keV^-2 to cm^2 


T_max=lambda En:(2*En*En)/(M+En+En)
C=lambda En,T_th:nat_unit_fac*quad(lambda T: C_T(En,T),T_th,T_max(En))[0]
#^ in units of cm^2 units

#Multiply by appropriate units to convert from natural units (keV^-2) to SI units (cm^2)

t=1#sec
#N=lambda En: t*n_flux(En)*quad(lambda T: norm_n_flux(En)*((Gf2*M/(2*math.pi))*(((G_V+G_A)**2)+(((G_V-G_A)**2)*((1-(T/En))**2))+(((G_V**2)-(G_A**2))*(M*T/(En*En))))),T_min,T_max(En))[0]
N=lambda En,T_th: t * total_flux * norm_n_flux(En) * (M_det/M) * C(En,T_th)

A=lambda T_th:quad(lambda En:N(En,T_th),E_min(T_th),np.inf)[0]
norm_N=lambda En,T_th:N(En,T_th)/A(T_th)
norm_N_th=lambda T_th: lambda En: norm_N(En,T_th)


E_min=lambda T_th: 0.5*(T_th+(((2*M*T_th)+(T_th*T_th))**0.5))

N_th=lambda T_th: quad(lambda En:N(En,T_th),E_min(T_th),np.inf)[0] * 1.4784374143281411 *60*60*24 # Correction to account for Joel's (1 day)

N_delE_th=lambda T_th,E1,E2: quad(lambda En:N(En,T_th),E1,E2)[0] * 1.4784374143281411 # Correction to account for Joel's
#th: Nuclear Recoil Energy(keV)

C_T1_T2=lambda T1,T2:quad(lambda T: quad(lambda E:nat_unit_fac*t*total_flux*norm_n_flux(E)*C_T(E,T)*(M_det/M)*1.4784374143281411,E_min(T),20000)[0],T1,T2)[0]
CS_T1_T2=lambda T1,T2:quad(lambda T: quad(lambda E:nat_unit_fac*C_T(E,T)*1.4784374143281411,E_min(T),20000)[0],T1,T2)[0]

import random



"""            
#E=np.arange(1e2,20e3,1e2)
#SENSITIVITY STUDY:-

import time

#draw_random_number_from_pdf(norm_N_th(0.01),[0,20e3])

def pdf_th(T_th):
    A_th=A(T_th)
    E=np.arange(E_min(T_th),20e3,50)
    p=[(E[1]-E[0])*N(e,T_th)/A_th for e in E]
    return E,p



def sim_data(events, threshold):   
    E,pdf=pdf_th(threshold)  
    pdf=np.fabs(pdf)
    pdf=pdf/sum(pdf)
    dist=np.random.choice(E,events,p=pdf) 
    return dist,np.array(E),np.array(pdf)
    
N_data=int(N_th(th))
data_10eV,E,pdf=sim_data(N_data, th)

print(len(E))
sim_dist=np.histogram(data_10eV, bins=np.insert(E,0,0))[0]
print(len(E))
P=lambda a,m2,E: a*(np.sin(m2*1.27*((L/1e3)/(E/1e6)))**2)#length from m to km, energy from keV to GeV
oscillated_pdf=lambda a,m2: pdf*(1-P(a,m2,E))
oscillated_df=np.vectorize(lambda a,m2: oscillated_pdf(a,m2)*N_data)

from scipy.stats import chisquare

#print(chisquare(sim_dist,f_exp=,ddof=len(E)))
start_time=time.time()
adelta=0.01
m2delta=0.01
a=np.arange(0+adelta, 1+adelta, adelta)#dimensionless
m2=np.arange(0+m2delta, 10+m2delta, m2delta)#eV^2
As,M2s=np.meshgrid(a, m2)
from matplotlib import colors, ticker, cm#for contourf ()
#v=lambda a,m2: sum(((oscillated_df(a,m2)-sim_dist)**2)/oscillated_df(a,m2))
v=lambda a,m2: chisquare(sim_dist,f_exp=oscillated_df(a,m2),ddof=len(E))[0]
chi_vec=np.vectorize(v)
V=chi_vec(As,M2s)
import scipy
cs=plt.contour(As,M2s,V, levels=scipy.stats.chi2.ppf(q = [0.9,0.95], df = len(E)), colors=['b','r'])#C.Ls corresponding to dof=389
lines=[cs.collections[0],cs.collections[1]]
labels = ['90% C.L.','95% C.L.']
plt.legend(lines, labels, loc=3)
plt.xscale("log")
plt.yscale("log")
plt.ylim(1e-1,10)
plt.xlim(1e-1,1)
plt.xlabel("sin^2(2*mixing angle)")
plt.ylabel("Mass Square Difference (eV^2)")
plt.title("Sensitivity @ "+str(L)+"m from reactor\nRun time="+str(t/tyr)+"years("+str(N_data)+" events)\nDetector threshold= "+str(th*1e3)+"eV")
print("--- %s seconds ---" % (time.time() - start_time)) 
"""


"""
start_time=time.time()
sim_data=lambda events, T_th: [E[weighted_choice_sub(pdf_th(T_th))] for i in range(events)]
data_1000_10eV=sim_data(10, 0.01)
print("--- %s seconds ---" % (time.time() - start_time)) 
"""    


"""
from scipy import stats
norm_N_e=lambda T_th: np.vectorize(norm_N)(E,T_th)
p=norm_N_e(0.01)
#rand_gen=stats.rv_discrete(name='rand_gen', values=(E, norm_N_e(0.01)))
#sim_data=rand_gen.rvs(size=100)
"""
"""
print("Recoil Analysis:-")
#----

"""

"""
#Smearing Matrix
def pad(S, pad_value=0):
    h=len(S[-1])    
    for i in range(len(S)):
        S[i]=S[i]+((h-len(S[i]))*[pad_value])
    return S

dist_f=lambda En,T: norm_n_flux(En) * C_T(En,T)
norm_fac= lambda En: quad(lambda T: dist_f(En,T), 0.01, T_max(En))[0]
dist=lambda En,T: dist_f(En,T)/(norm_fac(En))
dist_vec=np.vectorize(dist)

import time
print("Build Smearing Matrix")
start_time=time.time()
T_th=0.01

dt=0.01
de=10

Smat=[[dist_f(e,r) for r in np.arange(T_th,T_max(e),dt)] for e in np.arange(E_min(T_th),20e3,de)]
Smat=pad(Smat,0)
Smat=np.array(Smat)
Smat=[de*Smat[i]/(sum(Smat[i])*dt) for i in range(len(Smat))]
Smat=np.array(Smat)
Smat=np.nan_to_num(Smat)
print("--- %s seconds ---" % (time.time() - start_time))
print(Smat.shape)


#from matplotlib import ticker

#ratio_mean_dist=lambda En, T_th: (1/En)*quad(lambda T: dist(En,T)*T,T_th,T_max(En))[0]/quad(lambda T: dist(En,T),T_th,T_max(En))[0]
"""
"""
es=lambda T_th: np.arange(E_min(T_th),20e3,100)
l1,=plt.plot(es(0.01),np.vectorize(ratio_mean_dist)(es(0.01),0.01))
l2,=plt.plot(es(0.05),np.vectorize(ratio_mean_dist)(es(0.05),0.05))
l3,=plt.plot(es(0.1),np.vectorize(ratio_mean_dist)(es(0.1),0.1))
plt.legend([l1,l2,l3],["Threshold=10eV","Threshold=50eV","Threshold=100eV"], loc=4)
plt.xlabel("Incident Neutrino Energy")
plt.ylabel("Mean Recoil Energy/Neutrino Energy")
#plt.yscale("log")
"""

"""
r=np.arange(0.01,T_max(20e3),0.01)
x=np.arange(E_min(0.01),20e3,100)
X,R=np.meshgrid(x,r)
X_R=dist_vec(X,R)
CS=plt.contourf(X,R,X_R, locator=ticker.LogLocator());plt.colorbar(); plt.yscale("log")
#plt.clabel(CS, inline=1, fontsize=10)
plt.xlabel("Incident Neutrino Energy")
plt.ylabel("Recoil Energy (keV)")
plt.title("")
"""

"""
plt.xlabel("Recoil Energy (keV)");plt.ylabel("ratio b/w dN/dT & E");plt.yscale("log"); plt.xscale("log")

ts=lambda En: np.arange(0.01,T_max(En),0.01)
l1,=plt.plot(ts(5e3),dist(5e3,ts(5e3)))
l2,=plt.plot(ts(10e3),dist(10e3,ts(10e3)))
l3,=plt.plot(ts(15e3),dist(15e3,ts(15e3)))
l4,=plt.plot(ts(20e3),dist(20e3,ts(20e3)))

plt.legend([l1,l2,l3,l4],["E_nu=5MeV","E_nu=10MeV","E_nu=15MeV","E_nu=20MeV"])
"""
"""
R=np.arange(1e-2,(1)+1e-2,1e-2)
n=[N_th(r) for r in R]
plt.plot(R*1e3,n); plt.xscale("log"); 
plt.xlabel("Detector Threshold Energy (eV)")
plt.ylabel("# of Events/second")
plt.title("Number of Coherent Events detected vs. Detector Threshold Energy")
plt.yscale("log")
"""

#N_th=lambda T_th: trapz([N(e,T_th) for e in np.arange(E_min(T_th),20e3,1)], dx=1)#E=E_min - 20MeV with 1keV step 




"""
n=[]
for e in np.arange(1,20001,100):#keV
    print(e)    
    temp=[]    
    for r in np.arange(1e-3,(1)+1e-3,1e-3):#eV
        #print("\t",r)        
        temp.append(N(e,r))
    n.append(temp)

from matplotlib.colors import LogNorm

n=np.flipud(np.array(n).transpose())
plt.imshow(n,norm=LogNorm(),extent=[0,1,0,1])
plt.colorbar()

plt.xlabel("Energy of Anti-Neutrino")
plt.ylabel("Detector Threshold Energy")
"""
#plt.yscale("log")
#plt.xscale("log")


"""
E=np.arange(1,(20e3)+1,1)


print(total_n)
#print(quad(lambda En:N(En,T_th),1,20e3+1)[0])

flux=[n_flux(e) for e in E]
l,=plt.plot(E,flux)
plt.yscale("log")

E=np.arange(E_min(1e-2),(20e3)+1,1)
n1=[N(e,1e-2) for e in E]
l1,=plt.plot(E,n1)

E=np.arange(E_min(1e-1),(20e3)+1,1)
n2=[N(e,1e-1) for e in E]
l2,=plt.plot(E,n2)

plt.legend([l,l1,l2],['Anti-electron Neutrino Spectrum','Observed in Detector with threshold=10eV','Observed in Detector with threshold=100eV'])

plt.xlabel("Incident Neutrino Energy (keV)")
plt.ylabel("Number of Events per sec")
plt.title("Number of Coherent Neutrino-Nucleus Interactions")
"""