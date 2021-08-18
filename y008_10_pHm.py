'''
Created on Sep 21, 2020

@author: keiin
'''

from pylab import *
from FigSetting2 import *
from sf import *

pHmin = 2
pHmax = 10
pH = arange(pHmin,pHmax+1e-10,0.01)

def Cfixcal0():
    rho = 997.102/1000 #(kg*1000/m3)(kg/L) Water Density (http://www.csgnetwork.com/h2odenscalc.html)
    
    DIC = arange(0,2e3+1e-10,1e1) #(umol/L)
    DIC = 993   #(Burns 1987)
    for i in arange(size(pH)):  #This is to pick i value where pH = 7.59, which is used for reference based on Burns 1987
        if pH[i]>7.59-1e-10 and pH[i]<7.59+1e-10:
            ii = i
    
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
    

    Cfix1 = CO2/CO2[ii]
    
    Cfix2 = CO2/(CO2 + Kco2)
    Cfix2 = Cfix2/Cfix2[ii]
    
    return Cfix1,Cfix2,CO2

def Cfixcal(VmaxD):
    rho = 997.102/1000 #(kg*1000/m3)(kg/L) Water Density (http://www.csgnetwork.com/h2odenscalc.html)
    
    DIC = arange(0,2e3+1e-10,1e1) #(umol/L)
    DIC = 993   #(Burns 1987)
    for i in arange(size(pH)):
        if pH[i]>7.59-1e-10 and pH[i]<7.59+1e-10:
            ii = i
    
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
    
    return Cfix,CO2in


CfixMax,CfixMin,CO2 = Cfixcal0()
VmaxD = arange(0,1500+1e-10,300)
VmaxD[0] = 0.1
VmaxD = array([150, 300, 500, 700, 950])
#VmaxD = VmaxD[::-1]
Colors = ['#000096','#1D24D3','#3A49DB','#576DE2','#7592E9','#92B6F0','#AFDBF8']
Colors = Colors[::-1]

Cfix = zeros((len(VmaxD),len(pH)))*nan
CO2in = copy(Cfix)
print(shape(Cfix))

figure(100)
plot([7.59,7.59],[0,50],':',color='k')
plot(pH,CfixMin,label = 'V$_{max}$ << D',color=Colors[0])
for i in arange(len(VmaxD)):
    Cfix[i,:],a = Cfixcal(VmaxD[i])
    if VmaxD[i]<1:
        Label = str(VmaxD[i])
    else:
        Label = 'V$_{max}$/D='+str(int(VmaxD[i]))
    plot(pH,Cfix[i],label = Label,color=Colors[i+1])
plot(pH,CfixMax,label = 'V$_{max}$ >> D',color=Colors[-1])
#fill_between(pH, CfixMax, CfixMin, where = CfixMax>=CfixMin, facecolor='#BDD7EE',edgecolor = "none")
xlabel('pH$_{m}$')
ylabel('Relative C fixation rate')
xlim(pHmin,pHmax)
ylim(0,65)

current_handles, current_labels = gca().get_legend_handles_labels()
reversed_handles = list(reversed(current_handles))
reversed_labels = list(reversed(current_labels))
legend(reversed_handles,reversed_labels,edgecolor='k',fontsize=18,facecolor='white',framealpha=1)

sf('','various',300)


Colors = ['#680000','#7E1D16','#933A2C','#A95742','#BE7557','#D4926D','#E9AF83']
Colors = Colors[::-1]

figure(101)
#plot([7.59,7.59],[0,20],':',color='k')
plot(pH,CO2,label = 'V$_{max}$ << D',color=Colors[0])
for i in arange(len(VmaxD)):
    a,CO2in[i,:] = Cfixcal(VmaxD[i])
    if VmaxD[i]<1:
        Label = str(VmaxD[i])
    else:
        Label = 'V$_{max}$/D='+str(int(VmaxD[i]))
    plot(pH,CO2in[i],label = Label,color=Colors[i+1])
plot(pH,CO2*0,label = 'V$_{max}$ >> D',color=Colors[-1])

xlabel('pH$_{m}$')
ylabel('[CO$_{2}$]$_{p}$ ($\mu$mol L$^{-1}$)')
xlim(pHmin,pHmax)
ylim(-20,1050)

current_handles, current_labels = gca().get_legend_handles_labels()
reversed_handles = list(reversed(current_handles))
reversed_labels = list(reversed(current_labels))
legend(reversed_handles,reversed_labels,edgecolor='k',fontsize=18)

sf('','CO2 various',300)

show()





