# FYS5555 PROJECT 1

import numpy as np 
import pandas as pd 
import scipy.integrate as integrate 
import matplotlib.pyplot as plt

# functions 

def diff_CS(s,M):
	return 3* 1/(2.56810e-9)*2*np.pi/(64*np.pi**2*s)*p1/p * M

#--------------------------------------------------------------

# Import comphep data 
TotalCS_comphep = np.loadtxt("comphep/SM_CS.txt",skiprows=3)
Asymmetry_comphep = np.loadtxt("comphep/SM_Asymm.txt",skiprows=3)


# Defining constants
q = 1
Q = 0.3333
e = 0.31345

sinW = 0.48076
sin2W = sinW**2
thetaW = np.arcsin(sinW)

gZ = e/(sinW*np.cos(thetaW))
g = e/sinW

ca = 0.5*gZ*(-0.5)
cv = -0.5*gZ*0.04
Delta = cv**2-ca**2
Sigma = cv**2+ca**2



m = 0.1057 #GeV
m_square = m**2


mZ = 91.1876
WidthZ = 2.43631
mZ_square = mZ**2

mH = 125.09
WidthH = 0.0061744
mH_square = mH**2

mW = mZ * np.cos(thetaW)
mW_square = mW**2


ee = (0.31345)**2 #e^2
#---------------------------------------------------------------

# Setting up model

CMenergy = np.linspace(10, 200, 200)

TotalCS = np.zeros(np.size(CMenergy))
TotalCS_1 = np.zeros(np.size(CMenergy))
TotalCS_2 = np.zeros(np.size(CMenergy))
TotalCS_3 = np.zeros(np.size(CMenergy))

Fw = np.zeros(np.size(CMenergy))
Bw = np.zeros(np.size(CMenergy))
Fw_1 = np.zeros(np.size(CMenergy))
Bw_1 = np.zeros(np.size(CMenergy))
Fw_2 = np.zeros(np.size(CMenergy))
Bw_2 = np.zeros(np.size(CMenergy))

X = np.linspace(-1,1,500)


Masses  = [4.85]#[0.511, 1.29, 4.85] 


M = 4.85#Masses[j]

M_square = M**2

Ca1 = [-0.5, 0-5, -0.5]
Cv1 = [-0.04, 0.19, -0.35]
ca1 = 0.5*gZ*(-0.5)#Ca1[j]
cv1 = 0.5*gZ*(-0.35)#Cv1[j]
Delta1 = cv1**2-ca1**2
Sigma1 = cv1**2+ca1**2


# Loop over energies

for i in range(len(CMenergy)):

	# Kinematics

	Ecm = CMenergy[i] #100.0 #GeV
	s = Ecm**2
	E_square = (0.5*Ecm)**2
	p = np.sqrt( E_square - m_square )
	p1 = np.sqrt( E_square - M_square )


	pk = (E_square + p**2)*np.ones(np.size(X))
	#pp1 = E_square + p*p1*X
	pp1 = E_square - p*p1*X
	#pk1 = E_square - p*p1*X
	pk1 = E_square + p*p1*X
	#kp1 = E_square - p*p1*X
	kp1 = E_square + p*p1*X
	#kk1 = E_square + p*p1*X
	kk1 = E_square - p*p1*X
	p1k1  = (E_square + p1*p1)*np.ones(np.size(X))


	# Amplitudes

	M_AA = ee**2/(9*s**2) * 8*( kp1*pk1 + kk1*pp1 + M_square*pk + m_square*p1k1 + 2*M_square*m_square*np.ones(np.size(X)) )

	k_ZZ = 8*np.real( 1/( (s-mZ_square-1j*mZ*WidthZ)*(s-mZ_square+1j*mZ*WidthZ) ) )
	M_ZZ = k_ZZ*( Sigma*Sigma1*(kk1*pp1+pk1*kp1) + Delta1*Sigma*M_square*pk + Delta*Sigma1*m_square*p1k1 + 2*Delta*Delta1*M_square*m_square*np.ones(np.size(X))-4*ca*cv*ca1*cv1*( pp1*kk1-kp1*pk1 ))



	M_ZA = np.real(ee/(3*s*(s-mZ_square+1j*mZ*WidthZ))) * \
		8*( cv*cv1*(kp1*pk1 + pp1*kk1 + m_square*p1k1 + M_square*pk + 2*m_square*M_square*np.ones(np.size(X))) + ca*ca1*(pk1*kp1-pp1*kk1) )		


	M_HH = 2.9e-10*np.real(1/( (s-mH_square-1j*mH*WidthH)*(s-mH_square+1j*mH*WidthH) )) * (pk*p1k1 - M_square*pk - m_square*p1k1 + M_square*m_square*np.ones(np.size(X)))
	# Due to problems of numerical precision, the prefactor has been calculated separately


	M_ZH = np.real ( \
		(( (s-mH_square-1j*mH*WidthH)*(s-mZ_square+1j*mZ*WidthZ) + (s-mH_square+1j*mH*WidthH)*(s-mZ_square-1j*mZ*WidthZ) )*4*mW_square/(ee*g**2*m*M) )**(-1) ) *4*M*m*cv*cv1*\
		(kk1-kp1-pk1+pp1)

	M_AH = np.real( ( 12*mW_square*s*(s-mH_square+1j*WidthH*mH)/(ee*g**2*m*M) )**(-1) )*\
		4*M*m*(kk1-kp1-pk1+pp1)



	M = M_ZZ + M_HH + M_AA + 2*M_ZA + 2*M_AH + M_ZH

	# if i==0:
	# 	factor = 3
	# else:
	# 	factor = 1
	diffCS = diff_CS(s,M)#/factor
	diffCS_A = diff_CS(s,M_AA+M_ZZ)
	#diffCS_Z = diff_CS(s,M_ZZ)
	#diffCS_H = diff_CS(s,M_HH)


	# plotting the differential CS
	#plt.plot(X,diffCS)
	#plt.plot(X,diffCS_A,linestyle='--')
	#plt.plot(X,diffCS_H,linestyle='--')
	#plt.plot(X,diffCS_AZ,linestyle='--')
	#plt.ylim([0.8, 1.2])
	#plt.xlabel(r'$\cos\theta$')
	#plt.ylabel(r'd$\sigma$/d$\cos\theta$ [pb/rad]')
	#plt.grid()
	#plt.savefig('figures/150gev_termH')
	#plt.show()


	# Calculating the total cross section
	TotalCS[i] = integrate.simps(diffCS,X)
	#TotalCS_1[i] = integrate.simps(diffCS_A,X)
	#TotalCS_2[i] = integrate.simps(diffCS_Z,X)
	#TotalCS_3[i] = integrate.simps(diffCS_H,X)


	Bw[i] = integrate.simps(diffCS[0:250])
	Fw[i] = integrate.simps(diffCS[251:])
	Bw_1[i] = integrate.simps(diffCS_A[0:250])
	Fw_1[i] = integrate.simps(diffCS_A[251:])
	#Bw_2[i] = integrate.simps(diffCS_AH[0:250])
	#Fw_2[i] = integrate.simps(diffCS_AH[251:])


#plt.title(r'$\sqrt{s}=10$ GeV')
#axs= plt.gca()
#axs.legend([r'total',r'$M_{\gamma}^2$'])#+ M_Z^2 + 2M_{\gamma Z}$',r'$M_{H}^2 + 2M_{\gamma}M_H + M_{ZH}$'])
#plt.savefig('figures/10gev_diff_A')
#plt.show()

#----------------------------------------------------------------------	

# plot total cross section
#plt.figure()
#plt.semilogy(CMenergy, TotalCS, linestyle='none', marker='.')
#plt.semilogy(TotalCS_comphep[:,0],TotalCS_comphep[:,1], linewidth=0.5)
#plt.grid('on',linestyle='--')
#plt.title(r'Total cross section $\mu \bar{\mu} \rightarrow b \bar{b}$')
#plt.xlabel(r'$\sqrt{s}$ [GeV]')
#plt.ylabel(r'$\sigma$ [pb]')
#axs = plt.gca()
#axs.legend(['Analytical','CompHEP'])
#plt.savefig('figures/TotalCS_CH.jpg')
#plt.show()

# plot asymmetry factor
Asymm = (Fw - Bw)/(Fw + Bw)
Asymm_A = (Fw_1 - Bw_1)/(Fw_1 + Bw_1)
#Asymm_ZA = (Fw_2 - Bw_2)/(Fw_2 + Bw_2)

#plt.figure()
plt.plot(CMenergy, Asymm, linestyle='-')
plt.plot(CMenergy, Asymm_A, linestyle='-')
#plt.plot(Asymmetry_comphep[:,0], Asymmetry_comphep[:,1], linewidth=0.5)

plt.grid('on',linestyle='--')
plt.title(r'Asymmetry factor')
plt.xlabel(r'$\sqrt{s}$ [GeV]')
plt.ylabel(r'$F$')
axs = plt.gca()
#axs.legend(['Analytical','CompHEP'])
axs.legend([r'$M_{\gamma}^2$+ $M_Z^2$ + $2M_{\gamma Z}$',r'$M_{\gamma}^2$+ $M_Z^2$'])
plt.savefig('figures/Asymm_noint.jpg')
plt.show()


#------------------------------------------------------------


