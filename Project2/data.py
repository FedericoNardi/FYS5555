import numpy as np 
import matplotlib.pyplot as plt


# LEAKAGE TEST
# Inflating curve - get the state equation
p = np.array([0.73, 0.99, 1.20, 1.41, 1.60, 1.76]) #mbar
dp = 0.01*np.ones(np.size(p))
V = np.array([100, 200, 300, 400, 500, 600]) #mL

coeff = np.polyfit(p,V,2)

plt.figure()
Vpred = coeff[2]+p*coeff[1]+p**2*coeff[0]
plt.errorbar(p,V,xerr=dp,marker='.',linestyle='none',elinewidth=0.7,capsize=0.5)
plt.plot(p,Vpred,'--',linewidth=0.5)
plt.grid(linestyle='--')
plt.legend(['fit','data'])
plt.title(r'Leakage test - inflating curve')
plt.xlabel(r'p [mbar]')
plt.ylabel(r'V [mL]')
plt.savefig('figures/inflate.eps')
#plt.show()

# leakage test
p = np.array([1.79, 1.64, 1.52, 1.42, 1.32, 1.24, 1.12, 1.06, 0.97, 0.91, 0.86, 0.78, 0.74, 0.67, 0.65])
dp = 0.02*np.ones(np.size(p))
T = np.array([27, 27, 27.1, 27.1, 27.1, 27.2, 27.2, 27.2, 27.2, 27.3, 27.3, 27.4, 27.4, 27.4, 27.3])
Tinfl = 26.7
# predict volume in the chamber
V = coeff[0]*p**2 + coeff[1]*p + coeff[2]
dV = (coeff[0]*2*p + coeff[1])*dp
# correct volume for temperature changes
V = V*T/Tinfl
t = np.linspace(0,70,15)

plt.figure()
plt.errorbar(t,V,yerr=dV,linestyle='--',marker='o',linewidth=0.5,elinewidth=0.7,capsize=0.5)
plt.grid(linestyle='--')
plt.ylabel(r'V [mL]')
plt.xlabel(r't [min]')
plt.title(r'Chamber leakage')
plt.savefig('figures/leak.eps')
#plt.show()

# Select final points to get leakage
tleak = t[10:]
Vleak = V[10:]
# fit to get the slope
coeff_leak = np.polyfit(tleak, Vleak, 1)

plt.figure()
plt.errorbar(tleak, Vleak, yerr=dV[10:],marker='o', linestyle='none')
plt.plot(tleak, coeff_leak[0]*tleak+coeff_leak[1],'--',linewidth=0.5)
plt.legend(['fit - slope ='+str(round(0.06*coeff_leak[0],2))+r' [L/h]','data'])
plt.grid(linestyle='--')
plt.ylabel(r'V [mL]')
plt.xlabel(r't [min]')
plt.title(r'Chamber leakage - Linear fit')
plt.savefig('figures/leak_fit.eps')
#plt.show() 




# efficiency plateau - chamber 4
Vplus = np.array([5.92,6.25,6.5,6.75,6.97,7.24,7.5,7.76,8.01,8.25,8.5,8.76,8.99,9.24,9.49,9.75])
Vmin = np.array([6.02,6.26,6.52,6.75,7,7.25,7.49,7.76,8,8.25,8.49,8.75,9.01,9.25,9.48,9.75])
Voltage = Vplus+Vmin
Viplus = np.array([
	7.2,
	7.75,
	8.19,
	8.54,
	8.63,
	9.36,
	9.81,
	10.36,
	10.54,
	11.05,
	11.61,
	12.29,
	13.12,
	14.08,
	15.28,
	16.78
	])
Vimin = np.array([
	7.8,
	7.91,
	8.16,
	8.54,
	8.51,
	9.27,
	9.75,
	10.26,
	10.15,
	11,
	11.52,
	12.14,
	13,
	13.82,
	14.91,
	16.24
	])
Eff = np.array([
	0.012,
	0.038,
	0.105,
	0.196,
	0.265,
	0.412,
	0.503,
	0.559,
	0.612,
	0.647,
	0.670,
	0.692,
	0.706,
	0.716,
	0.735,
	0.729
	])
dEff = np.array([
	0.007,
	0.007,
	0.007,
	0.007,
	0.006,
	0.006,
	0.005,
	0.005,
	0.005,
	0.004,
	0.004,
	0.004,
	0.004,
	0.004,
	0.004,
	0.004
	])
iplus = Viplus - Vplus
imin = Vimin - Vmin

# Plot 
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel(r'V [kV]')
ax1.set_ylabel(r'$\eta$', color=color)
ax1.errorbar(Voltage, Eff,yerr=dEff, marker='o',linewidth=0.5, linestyle='--', markersize=4, color=color)
ax1.tick_params(axis='y', labelcolor=color)
plt.legend([r'efficiency $\eta$'],loc=2)
plt.grid(linestyle='--')
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel(r'$i$ [$\mu A$]', color=color)  # we already handled the x-label with ax1
ax2.set_ylim([0,20])
ax2.plot(Voltage, 0.5*(iplus+imin), marker='x',linewidth=0.5, linestyle='--', markersize=4, color=color)
plt.legend([r'dark current'],loc=4)
ax2.tick_params(axis='y', labelcolor=color)
plt.title('chamber test - efficiency plateau')
plt.savefig('figures/efficiency_plot.eps')
plt.show()