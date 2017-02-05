
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
#from CnNS import pdf_E
#from CnNS import random_1
#from CnNS import pdf_slice_1  
import sys as syst
from CNnsCS import *
#from CnNS_Counting_2 import *



#plt.plot(pdf_E,random_1)
#plt.plot(pdf_E,pdf_slice_1) 
fig, ax = plt.subplots()
# E in eV->keV

th=float(syst.argv[1])
s_b=float(syst.argv[2])
sys=float(syst.argv[3])
bins=float(syst.argv[4])
L=float(syst.argv[5])
m=float(syst.argv[6])

T=300#day
s=[C_T1_T2(i,i+((1-(th/1000))/bins))*(60*60*24)*T*m/(L*L) for i in np.linspace((th/1000), 1, num=bins+1)[:-1]]
pdf_E=np.linspace((th/1000), 1, num=bins+1)[:-1]
E_mid=(np.linspace((th/1000), 1, num=bins+1)[:-1]+np.linspace((th/1000), 1, num=bins+1)[1:])/2
pdf_slice_1=random_1=s

#print(len(pdf_E))
#print(len(random_1))
#print(len(pdf_slice_1))
  
def Oscillate(E_mid,s,a,m2):
    k=1.27 #Constant associated with Units   
    P=a*(np.sin(m2*k*((L/1e3)/(E_mid/1e6)))**2)#length from m to km, energy from keV to GeV
    pdf_1_osc=s*(1-P)#Energy Flux with Oscillations
    return pdf_1_osc

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
t = E_mid
a = 0.5
m2 = 1
plt.plot(E_mid*1e3,s,lw=2); 
#s = a0*np.sin(2*np.pi*f0*t)
s_o=Oscillate(E_mid,s,a,m2)
l ,=plt.plot(t*1e3, s_o, color='red', lw=1.5)
#plt.axis([0, 1, 0, 1])
plt.xlabel("Energy (eV)")
plt.ylabel("# of CNS Events")

#Postion of Slider
axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.28, 0.05, 0.65, 0.03], axisbg=axcolor)
axamp = plt.axes([0.28, 0.1, 0.65, 0.03], axisbg=axcolor)

sfreq = Slider(axfreq, r"sin$^2(2 \theta)$", 0, 1, valinit=a)
samp = Slider(axamp, r"$\Delta m^2 \ ($eV$^2)$", 1e-3, 10, valinit=m2)

def update(val):
    #m2 = np.log(samp.val)
    #samp.valtext.set_text(m2)
    m2 = samp.val
    a = sfreq.val
    l.set_ydata(Oscillate(pdf_E,s,a,m2))
    fig.canvas.draw_idle()
sfreq.on_changed(update)
samp.on_changed(update)

plt.title("Runtime  of "+str(T)+" days\nBlue: No Sterile\nRed: Yes Sterile")
plt.show()


