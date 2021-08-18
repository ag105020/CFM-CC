'''
Created on Sep 21, 2020

@author: keiin
'''

from pylab import *
from FigSetting2 import *
from sf import *

pH = arange(2,10.001,0.01)

def Cfixcal(VmaxD):
    rho = 997.102/1000 #(kg*1000/m3)(kg/L) Water Density (http://www.csgnetwork.com/h2odenscalc.html)
    
    DIC = arange(0,2e3+1e-10,1e1) #(umol/L)
    DIC = 993   #(Burns 1987)
    for i in arange(size(pH)):
        if pH[i]>7.59-1e-10 and pH[i]<7.59+1e-10:
            ii = i
    
    print(ii,pH[ii],pH[ii-1],pH[ii+1])
        
    
    H = 1/(10**pH) #(mol/L)
    H = H/rho #(mol/kg)
    
    T = 25 + 273.15
    S = 35 
    Kco2 = 44 #(uM) from Yoshizawa-san's email (Jensen 2020)
    
    pK1 = 3633.86/T - 61.2173 + 9.6777*log(T) - 0.011555*S + 0.0001152*S**2 #(Emerson book, pp. 131)
    pK2 = 471.78/T + 25.9290 - 3.16967*log(T) - 0.01781*S + 0.0001122*S**2  #(Emerson book, pp. 131)
    
    K1 = 1/(10**(pK1)) #(mol/kg)
    K2 = 1/(10**(pK2)) #(mol/kg)
    
    CO2 = DIC/(1 + K1/H + K1*K2/H**2) #(mol/L)
    HCO3 = DIC/(H/K1 + 1 + K2/H)  #(mol/L)
    CO3 = DIC/(1 + H**2/(K1*K2) + H/K2) #(mol/L)
    
    #===================
    # Inner CO2
    #===================
    #VmaxD = 500  #Vmax/D
    
    a = 1
    b = VmaxD + Kco2 - CO2
    c = -CO2*Kco2
    
    CO2in = (-b + (b**2 - 4*a*c)**0.5)/(2*a)
    #======================
    # Photosynthesis rate
    #======================
    Cfix = CO2-CO2in
    Cfix = Cfix/Cfix[ii]
    #======================
    Cfix1 = CO2/CO2[ii]
    
    Cfix2 = CO2/(CO2 + Kco2)
    Cfix2 = Cfix2/Cfix2[ii]
    
    return Cfix,Cfix1,Cfix2,CO2

VmaxD = 280
Cfix,Cfix1,Cfix2,CO2 = Cfixcal(VmaxD) 

figure(0)
plot(pH,Cfix)
xlabel('pH$_{m}$')
ylabel('C fixation rate (dimensionless)')
title('V$_{max}$/D = '+str(VmaxD),y=1.02)

figure(1)
plot(pH,Cfix1)
xlabel('pH$_{m}$')
ylabel('C fixation rate (dimensionless)')
title('V$_{max}$ >> D', y=1.02)
   
figure(2)
plot(pH,Cfix2)
xlabel('pH$_{m}$')
ylabel('C fixation rate (dimensionless)')
title('V$_{max}$ << D', y=1.02)


# 
# figure(3)
# plot(pH,Cfix,label='DIC='+str(DIC))
# xlabel('pH')
# ylabel('C fixation rate (dimensionless)')
# title('Full solution', y=1.02)
# 
figure(4)
plot(pH,CO2,color='#680000')
xlabel('pH$_{m}$')
ylabel('[CO$_{2}]_{m}$ ($\mu$mol L$^{-1}$)')
xlim(2,10)
ylim(-20,1050)
sf('','CO2m',300)
#title('[CO$_{2}$]', y=1.02)
# 
# figure(5)
# plot(pH,CO2in,label='DIC='+str(DIC))
# xlabel('pH')
# ylabel('CO$_{2}$ ($\mu$M)')
# ylim(top=1033)
# title('[CO$_{2}$]$_{in}$', y=1.02)

show()





