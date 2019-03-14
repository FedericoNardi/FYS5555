import numpy as np 
import pandas as pd 
import scipy.integrate as integrate 
import matplotlib.pyplot as plt

# Import comphep data 
# Ecm = [1000, 2500]#gev

# for i in range(len(Ecm)):
# 	diffCS = np.loadtxt('comphep/Zprime_diffCS_'+str(Ecm[i])+'gev.txt',skiprows=3)
# 	plt.plot(diffCS[:,0],diffCS[:,1])
# plt.grid()
# axs = plt.gca()
# axs.legend([r'$\sqrt{s} = $'+str(Ecm[0]/1e3)+' TeV',r'$\sqrt{s} = $'+str(Ecm[1]/1e3)+' TeV'])
# plt.xlabel(r'$\cos\theta$')
# plt.ylabel(r'd$\sigma$/d$\cos\theta$ [pb/rad]')
# plt.savefig('figures/Zprime_diffCS')
# plt.show()

# TotalCS = np.loadtxt('comphep/Zprime_CS.txt',skiprows=3)
# plt.figure()
# plt.semilogy(TotalCS[:,0],TotalCS[:,1])
# plt.grid()
# plt.xlabel(r'$\sqrt{s}$ [GeV]')
# plt.ylabel(r'$\sigma$ [pb]')
# plt.savefig('figures/Zprime_TotalCS')
# plt.show()


# Asymm = np.loadtxt('comphep/Zprime_asymm.txt',skiprows=3)
# plt.figure()
# plt.plot(Asymm[:,0],Asymm[:,1])
# plt.grid()
# plt.xlabel(r'$\sqrt{s}$ [GeV]')
# plt.ylabel(r'$F$')
# plt.savefig('figures/Zprime_asymm')
# plt.show()

################################################à
CrossSection_e = np.loadtxt('comphep/totalCS_e.txt',skiprows=3)
CrossSection_c = np.loadtxt('comphep/SM_CS.txt',skiprows=3)
CrossSection_b = np.loadtxt('comphep/totalCS_c.txt',skiprows=3)

Asymm_e = np.loadtxt('comphep/asymm_e.txt',skiprows=3)
Asymm_c = np.loadtxt('comphep/SM_Asymm.txt',skiprows=3)
Asymm_b = np.loadtxt('comphep/asymm_c.txt',skiprows=3)

plt.figure()
plt.semilogy(CrossSection_e[:,0],CrossSection_e[:,1])
plt.semilogy(CrossSection_c[:,0],CrossSection_c[:,1])
plt.semilogy(CrossSection_b[:,0],CrossSection_b[:,1])
plt.grid('on',linestyle='--')
plt.ylabel(r'$\sigma$ [pb]')
plt.xlabel(r'$\sqrt{s}$ [GeV]')
axs = plt.gca()
axs.legend([r'$\mu^+\mu^-\rightarrow e^+e^-$',r'$\mu^+\mu^-\rightarrow b \bar{b}$',r'$\mu^+\mu^-\rightarrow c \bar{c}$'])
plt.savefig('figures/Compare_CS')

plt.figure()
plt.plot(Asymm_e[:,0], Asymm_e[:,1])
plt.plot(Asymm_c[:,0], Asymm_c[:,1])
plt.plot(Asymm_b[:,0], Asymm_b[:,1])
plt.grid('on',linestyle='--')
plt.xlabel(r'$\sqrt{s}$ [GeV]')
axs = plt.gca()
axs.legend([r'$\mu^+\mu^-\rightarrow e^+e^-$',r'$\mu^+\mu^-\rightarrow b \bar{b}$',r'$\mu^+\mu^-\rightarrow c \bar{c}$'])
plt.savefig('figures/Compare_asymm')

plt.show()