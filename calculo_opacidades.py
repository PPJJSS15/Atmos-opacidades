# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 13:03:11 2022

"""
"Authors:"
"Pablo Jiménez Sánchez and Patxi Eguiguren Arrizabalaga"
import matplotlib.pyplot as plt
import numpy as np
#atmosferas
#%%
# EEE

def EcBoltzmann(N,g,U,chi,k,T):
    return N*(g/U)*np.exp(-chi/(k*T))

def EcSaha(ne,T,U1,U2,chi,k):
    return ne*(2.07e-16) * (U1/U2)* T**(-1.5) * np.exp(chi/(k*T))
#%% 
# Constantes

Teff = 5000
c = 3e8 #m/s

K = 8.617e-5 #eV/K
Kerg = K*1.60224e-12
h = 4.13566733e-15 #eV*s
R = 1.0968e-3 #Angstrom-1


Pe =8.5927
T = 5440.9
Pe =4.6795e+02
T = 8686.1
iter = 0
for T in [5440.9,8686.1,]:
    if T < 6000:
        Pe = 8.5927
    else :
        Pe = 4.6795e2
    # Propiedades especies
    g_HI=[2,8,18,32]
    Eexcit_HI=[0,10.199810,12.087494,12.748532]
    chi_HI=13.598434
    U_HI=2
    U_HII=1
    U_Hm= 1
    # chi_Hm = 14.359888
    chi_Hm = 0.755
    
    #%%
    # =============================================================================
    # # Poblaciones
    # =============================================================================
    # 1 Electrones
    ne = Pe/(Kerg*T)
    
    # 2 Saha
    NHI_NHII = EcSaha(ne,T,U_HI,U_HII,chi_HI,K)
    NHm_NHI = EcSaha(ne,T,U_Hm,U_HI,chi_Hm,K)
    
    # 3 Conservación Carga
    
    NHII = ne
    
    NHI = NHII * NHI_NHII
    NHm = NHI * NHm_NHI
    
    # A=[[0,1,-1],[1,-NHI_NHII,0],[-NHm_NHI,0,1]]
    # B = [ne,0,0]
    # solve = np.linalg.solve(A,B)
    # NHII =solve[1] 
    # NHI = solve[0] 
    # NHm =solve[2] 
    # 4 Boltzmann
    niveles_HI =[]
    
    for i in np.arange(0, 4):
        niveles_HI.append(EcBoltzmann(NHI , g_HI[i] , U_HI ,Eexcit_HI[i] , K , T ))
    print('NHII:  '+ str(NHII))
    print('NHI:  '+ str(NHI))
    print('NHm:  '+ str(NHm))
    print('ne:  '+ str(ne))
    
    print('---------------------------------')
    print('HI n=1:  '+str(niveles_HI[0]))
    print('HI n=2:  '+str(niveles_HI[1]))
    print('HI n=3:  '+str(niveles_HI[2]))
    print('HI n=4:  '+str(niveles_HI[3]))
    #%% 
    # =============================================================================
    # # Opacidades
    # =============================================================================
    lambdaa = np.arange(500,20000,10) #angstrom
    nu = c/(lambdaa*1e-10)
    lambda1m = 911.65
    lambda1p = 911.83 
    lambda2m = 3646.61 
    lambda2p =  3647.33
    lambda3m = 8204.86 
    lambda3p = 8206.51
    lambda4m = 14586.43 
    lambda4p = 14589.35
    lambdaHmm = 14999
    lambdaHmp = 15002
    #lambdaa = [lambda1m,lambda1p,lambda2m,lambda2p,lambda3m,lambda3p,lambda4m,lambda4p,lambdaHmm,lambdaHmp]
    #lambdaa = np.array([lambda1m,lambda1p,lambda2m,lambda2p,lambda3m,lambda3p,lambda4m,lambda4p,lambdaHmm,lambdaHmp])
    nu = c/(lambdaa*1e-10)
    # 1 Libre-Libre HI
    gaunt=1 + (0.3456/(lambdaa*R)**(1/3)) * ((lambdaa*1e-10*K*T/(h*c)) + 0.5)
    sigmalibrelibre_HI = (3.7e8)*(1./(nu**3)) * (1./(np.sqrt(T)))*gaunt
    opaclibrelibre_HI = sigmalibrelibre_HI*ne*NHII*(1.-np.exp(-(h*nu)/(K*T)))
    
    
    # 2 Libre-Libre H-
    f0 = -2.2763 -1.6850 * np.log10(lambdaa) + 0.76661 * np.log10(lambdaa)**2 - 0.053346* np.log10(lambdaa)**3
    f1 = 15.2827 - 9.1846* np.log10(lambdaa) + 1.99381 * np.log10(lambdaa)**2 - 0.142631* np.log10(lambdaa)**3
    f2 = -197.789 + 190.266* np.log10(lambdaa) - 67.9775* np.log10(lambdaa)**2 + 10.6913* np.log10(lambdaa)**3 - 0.625151* np.log10(lambdaa)**4
    theta = 5040/T
    sigmalibrelibre_Hm = 10**(-26) * 10**(f0+f1*np.log10(theta)+f2*np.log10(theta)**2)
    opaclibrelibre_Hm = Pe * sigmalibrelibre_Hm * NHI
    
    # 3 Ligado-Libre HI
    lambda_HI1 = 1/R
    lambda_HI2 = 4/R
    lambda_HI3 = 9/R
    lambda_HI4 = 16/R
    
    # 3.1 HI1
    counter = 0
    opacliglib_HI1 = np.zeros_like(opaclibrelibre_HI)
    opacliglib_HI2 = np.zeros_like(opaclibrelibre_HI)
    opacliglib_HI3 = np.zeros_like(opaclibrelibre_HI)
    opacliglib_HI4 = np.zeros_like(opaclibrelibre_HI)
    opacliglib_Hm = np.zeros_like(opaclibrelibre_HI)
    sigma_liglib_HI1 = np.zeros_like(opaclibrelibre_HI)
    sigma_liglib_HI2 = np.zeros_like(opaclibrelibre_HI)
    sigma_liglib_HI3 = np.zeros_like(opaclibrelibre_HI)
    sigma_liglib_HI4 = np.zeros_like(opaclibrelibre_HI)
    
    for i in np.arange(500,lambda_HI1 + 11,10):
        
        nu1 = c/(i*1e-10)
        gaunt=1 - (0.3456/(i*R)**(1/3)) * ((i*R/(1)) - 0.5)
        sigma_liglib_HI1[counter]=(2.815e29)*1/((nu1)**3)*gaunt
        opacliglib_HI1[counter] =sigma_liglib_HI1[counter]*niveles_HI[0]*((1)-np.exp(-h*nu1/(K*T)))
        
        counter+=1
        
    # 3.2 HI2
    counter = 0
    for i in np.arange(500,lambda_HI2 + 11,10):
        nu2 = c/(i*1e-10)
        gaunt=1 - (0.3456/(i*R)**(1/3)) * ((i*R/(4)) - 0.5)
        sigma_liglib_HI2[counter]=(2.815e29)*1/(32*(nu2)**3)*gaunt
        opacliglib_HI2[counter] =sigma_liglib_HI2[counter]*niveles_HI[1]*((1)-np.exp(-h*nu2/(K*T)))
        
        counter+=1
    # 3.3 HI3

    counter = 0
    for i in np.arange(500,lambda_HI3 + 11,10):
        nu3 = c/(i*1e-10)
        
        gaunt= 1 - (0.3456/(i*R)**(1/3)) * ((i*R/(9)) - 0.5)
        sigma_liglib_HI3[counter]=(2.815e29)*1/(243*(nu3)**3)*gaunt
        opacliglib_HI3[counter] =sigma_liglib_HI3[counter]*niveles_HI[2]*((1)-np.exp(-h*nu3/(K*T)))
        
        counter+=1
    
    # 3.4 HI4
    counter = 0
    for i in np.arange(500,lambda_HI4 + 11,10):
        nu2 = c/(i*1e-10)
        gaunt=1 - (0.3456/(i*R)**(1/3)) * ((i*R/(16)) - 0.5)
        sigma_liglib_HI4[counter]=(2.815e29)*1/(1024*(nu2)**3)*gaunt
        opacliglib_HI4[counter] =sigma_liglib_HI4[counter]*niveles_HI[3]*((1)-np.exp(-h*nu2/(K*T)))
        counter += 1
        
    figura = plt.figure(figsize=(7,5), dpi=200)
    plt.plot(lambdaa,sigma_liglib_HI1,label='n = 1')
    plt.plot(lambdaa,sigma_liglib_HI2,label='n = 2')
    plt.plot(lambdaa,sigma_liglib_HI3,label='n = 3')
    plt.plot(lambdaa,sigma_liglib_HI4,label='n = 4')
    plt.ylabel(r'$\sigma_{bf}  (cm^{2})$',fontsize = 13)
    plt.xlabel(r'$ \lambda (\AA) $',fontsize = 13)
    
    plt.xlim(0,15000)
    plt.legend()
    plt.show()
    # 4 Ligado-Libre Hm
    lambda_Hm = 15000
    a0 = 1.99654
    a1 = -1.18267e-5
    a2 = 2.64243e-6
    a3 = -4.40524e-10
    a4 = 3.23992e-14
    a5 = -1.39568e-18
    a6 = 2.78701e-23
    counter = 0
    for i in np.arange(500,lambda_Hm + 10,10):
        lambdaa = i
        sigma_liglib_Hm = (a0 + a1*lambdaa + a2*(lambdaa**2)+ a3*(lambdaa**3) + a4*(lambdaa**4) + a5*(lambdaa**5) + a6*(lambdaa**6))*1e-18
        opacliglib_Hm[counter] = 4.145e-10 * sigma_liglib_Hm * Pe * theta **(5/2) * 10**(0.754*theta) * NHI
        counter+=1
    # 5 Ligado-Ligado
    lambdaa = np.arange(500,20000,10)
    u=3
    l=2
    gbb = 0.869 - (3/2**3)
    f = 2**5 / (np.pi * 3**(1.5)) * (gbb/(l**5 * u**3)) * (1/l**2 - 1/u**2)**(-3)
    e = 4.8e-10
    m = 9.1094e10-28
    c = 3e8 * 1e2
    sigmabb = np.pi * f* e**2 /(m*c)
    # 5.1 Balmer alfa
    u=3
    l=2
    f = 0.6411
    sigmabb = np.pi * f* e**2 /(m*c)
    opacbb_Halpha = sigmabb * (niveles_HI[l-1]-niveles_HI[u-1]*g_HI[l-1]/g_HI[u-1])
    iBalfa = np.argwhere(lambdaa < 6563)[-1]
    # 5.2 Balmer beta
    u=4
    l=2
    f = 0.1194
    sigmabb = np.pi * f* e**2 /(m*c)
    opacbb_Hbeta = sigmabb * (niveles_HI[l-1]-niveles_HI[u-1]*g_HI[l-1]/g_HI[u-1])
    iBbeta = np.argwhere(lambdaa < 4861)[-1]
    # 5.3 Lyman alfa
    u=2
    l=1
    f = 0.4164
    sigmabb = np.pi * f* e**2 /(m*c)
    opacbb_Lalpha = sigmabb * (niveles_HI[l-1]-niveles_HI[u-1]*g_HI[l-1]/g_HI[u-1])
    iLalfa = np.argwhere(lambdaa < 1215)[-1]
    # 5.4 Lyman beta
    u=3
    l=1
    f = 0.6411
    sigmabb = np.pi * f* e**2 /(m*c)
    opacbb_Lbeta = sigmabb * (niveles_HI[l-1]-niveles_HI[u-1]*g_HI[l-1]/g_HI[u-1])
    iLbeta = np.argwhere(lambdaa < 1025)[-1]
    # 5.5 Lyman gamma
    u=4
    l=1
    f = 0.02901
    sigmabb = np.pi * f* e**2 /(m*c)
    opacbb_Lgamma = sigmabb * (niveles_HI[l-1]-niveles_HI[u-1]*g_HI[l-1]/g_HI[u-1])
    iLgamma = np.argwhere(lambdaa < 972)[-1] 
    # 5.6 Paschen alfa
    u=4
    l=3
    f = 0.8426
    sigmabb = np.pi * f* e**2 /(m*c)
    opacbb_Palpha = sigmabb * (niveles_HI[l-1]-niveles_HI[u-1]*g_HI[l-1]/g_HI[u-1])
    iPalfa = np.argwhere(lambdaa < 18750)[-1]
    # 6 Opacidad Electrones
    opac_el = 6.25e-25 * ne
    if iter <1:
        opacidad_total_5000 = opaclibrelibre_HI + opaclibrelibre_Hm + opacliglib_HI1 +opacliglib_HI2 +opacliglib_HI3 + opacliglib_HI4 + opacliglib_Hm +opac_el
        opacidad_total_5000[iBalfa] = opacidad_total_5000[iBalfa]+ opacbb_Halpha
        opacidad_total_5000[iBbeta] = opacidad_total_5000[iBbeta]+ opacbb_Hbeta
        opacidad_total_5000[iLalfa] = opacidad_total_5000[iLalfa]+ opacbb_Lalpha
        opacidad_total_5000[iLbeta] = opacidad_total_5000[iLbeta]+ opacbb_Lbeta
        opacidad_total_5000[iLgamma] = opacidad_total_5000[iLgamma]+ opacbb_Lgamma
        opacidad_total_5000[iPalfa] = opacidad_total_5000[iPalfa]+ opacbb_Palpha
    else:
        opacidad_total_8000 = opaclibrelibre_HI + opaclibrelibre_Hm + opacliglib_HI1 +opacliglib_HI2 +opacliglib_HI3 + opacliglib_HI4 + opacliglib_Hm +opac_el
        opacidad_total_8000[iBalfa] = opacidad_total_8000[iBalfa]+ opacbb_Halpha
        opacidad_total_8000[iBbeta] = opacidad_total_8000[iBbeta]+ opacbb_Hbeta
        opacidad_total_8000[iLalfa] = opacidad_total_8000[iLalfa]+ opacbb_Lalpha
        opacidad_total_8000[iLbeta] = opacidad_total_8000[iLbeta]+ opacbb_Lbeta
        opacidad_total_8000[iLgamma] = opacidad_total_8000[iLgamma]+ opacbb_Lgamma
        opacidad_total_8000[iPalfa] = opacidad_total_8000[iPalfa]+ opacbb_Palpha
    iter+=1
#%%
# =============================================================================
# GRÁFICOS
# =============================================================================

# 1 tau = 1 
lambdaa = np.arange(500,20000,10)
opacidad_total = opaclibrelibre_HI + opaclibrelibre_Hm + opacliglib_HI1 +opacliglib_HI2 +opacliglib_HI3 + opacliglib_HI4 + opacliglib_Hm +opac_el

figura = plt.figure(figsize=(7,5), dpi=200)
plt.plot(lambdaa,opacidad_total_5000,label='Teff = 5000 K')
plt.plot(lambdaa,opacidad_total_8000,label='Teff = 8000 K')
plt.ylabel(r'$\kappa_{total}  (cm^{-1})$',fontsize = 13)
plt.xlabel(r'$ \lambda (\AA) $',fontsize = 13)
plt.yscale('log')
plt.xscale('log')
#plt.ylim(0,1.1e-7)
plt.legend()
plt.show()

# =============================================================================
#  2º GRÁFICO
# =============================================================================
iter = 0
lambdaa = 3700
Temperatures = np.loadtxt('t5000.dat')[:,3]
Pes = np.loadtxt('t5000.dat')[:,4]
opacidad_total = np.zeros_like(Temperatures)
for T in Temperatures:
    Pe = Pes[iter]
    # Propiedades especies
    g_HI=[2,8,18,32]
    Eexcit_HI=[0,10.199810,12.087494,12.748532]
    chi_HI=13.598434
    U_HI=2
    U_HII=1
    U_Hm= 1
    # chi_Hm = 14.359888
    chi_Hm = 0.755
    
    #%%
    # =============================================================================
    # # Poblaciones
    # =============================================================================
    # 1 Electrones
    ne = Pe/(Kerg*T)
    
    # 2 Saha
    NHI_NHII = EcSaha(ne,T,U_HI,U_HII,chi_HI,K)
    NHm_NHI = EcSaha(ne,T,U_Hm,U_HI,chi_Hm,K)
    
    # 3 Conservación Carga
    
    NHII = ne
    
    NHI = NHII * NHI_NHII
    NHm = NHI * NHm_NHI
    
    # A=[[0,1,-1],[1,-NHI_NHII,0],[-NHm_NHI,0,1]]
    # B = [ne,0,0]
    # solve = np.linalg.solve(A,B)
    # NHII =solve[1] 
    # NHI = solve[0] 
    # NHm =solve[2] 
    # 4 Boltzmann
    niveles_HI =[]
    
    for i in np.arange(0, 4):
        niveles_HI.append(EcBoltzmann(NHI , g_HI[i] , U_HI ,Eexcit_HI[i] , K , T ))
    
    #%% 
    # =============================================================================
    # # Opacidades
    # =============================================================================
    #lambdaa = np.arange(500,20000,10) #angstrom
    nu = c/(lambdaa*1e-10)
    lambda1m = 911.65
    lambda1p = 911.83 
    lambda2m = 3646.61 
    lambda2p =  3647.33
    lambda3m = 8204.86 
    lambda3p = 8206.51
    lambda4m = 14586.43 
    lambda4p = 14589.35
    lambdaHmm = 14999
    lambdaHmp = 15002
    #lambdaa = [lambda1m,lambda1p,lambda2m,lambda2p,lambda3m,lambda3p,lambda4m,lambda4p,lambdaHmm,lambdaHmp]
    #lambdaa = np.array([lambda1m,lambda1p,lambda2m,lambda2p,lambda3m,lambda3p,lambda4m,lambda4p,lambdaHmm,lambdaHmp])
    
    # 1 Libre-Libre HI
    gaunt=1 + (0.3456/(lambdaa*R)**(1/3)) * ((lambdaa*1e-10*K*T/(h*c)) + 0.5)
    sigmalibrelibre_HI = (3.7e8)*(1./(nu**3)) * (1./(np.sqrt(T)))*gaunt
    opaclibrelibre_HI = sigmalibrelibre_HI*ne*NHII*(1.-np.exp(-(h*nu)/(K*T)))
    
    
    # 2 Libre-Libre H-
    f0 = -2.2763 -1.6850 * np.log10(lambdaa) + 0.76661 * np.log10(lambdaa)**2 - 0.053346* np.log10(lambdaa)**3
    f1 = 15.2827 - 9.1846* np.log10(lambdaa) + 1.99381 * np.log10(lambdaa)**2 - 0.142631* np.log10(lambdaa)**3
    f2 = -197.789 + 190.266* np.log10(lambdaa) - 67.9775* np.log10(lambdaa)**2 + 10.6913* np.log10(lambdaa)**3 - 0.625151* np.log10(lambdaa)**4
    theta = 5040/T
    sigmalibrelibre_Hm = 10**(-26) * 10**(f0+f1*np.log10(theta)+f2*np.log10(theta)**2)
    opaclibrelibre_Hm = Pe * sigmalibrelibre_Hm * NHI
    
    # 3 Ligado-Libre HI
    lambda_HI1 = 1/R
    lambda_HI2 = 4/R
    lambda_HI3 = 9/R
    lambda_HI4 = 16/R
       
    # 3.2 HI2
    nu2 = c/(lambdaa*1e-10)
    gaunt=1 - (0.3456/(lambdaa*R)**(1/3)) * ((lambdaa*R/(4)) - 0.5)
    sigma_liglib_HI2=(2.815e29)*1/(32*(nu2)**3)*gaunt
    opacliglib_HI2 =sigma_liglib_HI2*niveles_HI[1]*((1)-np.exp(-h*nu2/(K*T)))
        
    # 3.3 HI3
    nu3 = c/(lambdaa*1e-10)
        
    gaunt= 1 - (0.3456/(lambdaa*R)**(1/3)) * ((lambdaa*R/(9)) - 0.5)
    sigma_liglib_HI3=(2.815e29)*1/(243*(nu2)**3)*gaunt
    opacliglib_HI3 =sigma_liglib_HI3*niveles_HI[2]*((1)-np.exp(-h*nu3/(K*T)))

    # 3.4 HI4
    nu2 = c/(lambdaa*1e-10)
    gaunt=1 - (0.3456/(lambdaa*R)**(1/3)) * ((lambdaa*R/(16)) - 0.5)
    sigma_liglib_HI4=(2.815e29)*1/(1024*(nu2)**3)*gaunt
    opacliglib_HI4 =sigma_liglib_HI4*niveles_HI[3]*((1)-np.exp(-h*nu2/(K*T)))
     
    # 4 Ligado-Libre Hm
    lambda_Hm = 15000
    a0 = 1.99654
    a1 = -1.18267e-5
    a2 = 2.64243e-6
    a3 = -4.40524e-10
    a4 = 3.23992e-14
    a5 = -1.39568e-18
    a6 = 2.78701e-23
    sigma_liglib_Hm = (a0 + a1*lambdaa + a2*(lambdaa**2)+ a3*(lambdaa**3) + a4*(lambdaa**4) + a5*(lambdaa**5) + a6*(lambdaa**6))*1e-18
    opacliglib_Hm = 4.145e-10 * sigma_liglib_Hm * Pe * theta **(5/2) * 10**(0.754*theta) * NHI
    
    # 5 Ligado-Ligado
    
    u=3
    l=2
    gbb = 0.869 - (3/2**3)
    f = 2**5 / (np.pi * 3**(1.5)) * (gbb/(l**5 * u**3)) * (1/l**2 - 1/u**2)**(-3)
    e = 4.8e-10
    m = 9.1094e10-28
    c = 3e8 * 1e2
    sigmabb = np.pi * f* e**2 /(m*c)
    # 5.1 Balmer alfa
    u=3
    l=2
    f = 0.6411
    sigmabb = np.pi * f* e**2 /(m*c)
    opacbb_Halpha = sigmabb * (niveles_HI[l-1]-niveles_HI[u-1]*g_HI[l-1]/g_HI[u-1])
    
    # 6 Opacidad Electrones
    opac_el = 6.25e-25 * ne
    
    opacidad_total[iter] = opaclibrelibre_HI + opaclibrelibre_Hm  +opacliglib_HI3 + opacliglib_HI4 + opacliglib_Hm +opac_el 
    iter+=1
#%%
tau = np.loadtxt('t5000.dat')[:,1]
np.savetxt('opacidadtau3700.dat',opacidad_total)
figura = plt.figure(figsize=(7,5), dpi=200)
opacidad_total3600 = np.loadtxt('opacidadtau3600.dat')
opacidad_total3700 = np.loadtxt('opacidadtau3700.dat')
opacidad_total6562 = np.loadtxt('opacidadtau6562.dat')
plt.plot(tau,opacidad_total3600,label=r'$\lambda = 3600 \AA$')
plt.plot(tau,opacidad_total3700,label=r'$\lambda = 3700 \AA$',linestyle ='dashed',color='k')
plt.plot(tau,opacidad_total6562,label=r'$\lambda = 6562 \AA$')
plt.ylabel(r'$\kappa_{total}  (cm^{-1})$',fontsize = 13)
plt.xlabel(r'$ log(\tau_R)  $',fontsize = 13)
plt.yscale('log')

#plt.ylim(0,1.1e-7)
plt.legend()
plt.show()


#Pe =8.5927
#T = 5440.9
