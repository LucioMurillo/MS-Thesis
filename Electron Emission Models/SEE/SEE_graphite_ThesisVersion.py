import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spopt

def delta_ts(E,t1,t2,t3,t4,s):
    deltahat = deltahat_ts*(1 + t1*(1 - mu**t2))
    Ehat = Ehat_ts*(1 + t3*(1 - mu**t4))
    x = E/Ehat
    D = (s*x)/(s - 1 + x**s)
    return deltahat*D

# FP elastic yield + rediffused
def delta_e(Ep,Ehat_e,P1_inf,P1,W,p):
    exp = np.exp(-((np.abs(Ep - Ehat_e)/W)**p)/p)
    return P1_inf + (P1 - P1_inf)*exp

def delta_r(Ep,Er,P1_inf,r):
    exp = np.exp(-(Ep/Er)**r)
    return P1_inf*(1 - exp)

def delta_re(Ep,Ehat_e,P1e_inf,P1,W,p,Er,P1r_inf,r):
    e = delta_e(Ep,Ehat_e,P1e_inf,P1,W,p)
    r = delta_r(Ep,Er,P1r_inf,r)
    return e+r

# Graphite
Phi = 4.62 # work function [eV]

Farhang_data = np.loadtxt('../data/Graphite/Farhang_1993.txt', dtype=(str), delimiter=',').astype(float)
Pedgley_data = np.loadtxt('../data/Graphite/Pedgeley_1992.txt', dtype=(str), delimiter=',').astype(float)
Farhang_elastic_data = np.loadtxt('../data/Graphite/Farhang_1993_elastic.txt', dtype=(str), delimiter=',').astype(float)
Pedgley_elastic_data = np.loadtxt('../data/Graphite/Pedgley_1992_elastic.txt', dtype=(str), delimiter=',').astype(float) # Includes rediffused electrons

mu = 1
Ep = np.linspace(0,2000,100000)

deltahat_ts = np.max(Farhang_data[:,1])
Ehat_ts = Farhang_data[np.where(Farhang_data[:,1] == np.max(Farhang_data[:,1]))[0][0],0]
p_opt = 0.66,0.8,0.7,1,1.78
# delta_ts_opt,_ = spopt.curve_fit(delta_ts, Farhang_data[:,0], Farhang_data[:,1], p_opt,bounds=(0,[5,5,5,5,5]))
# delta_ts_poco = delta_ts(Ep,*delta_ts_opt)
delta_ts_poco = delta_ts(Ep, *p_opt)

deltahat_ts = np.max(Pedgley_data[:,1])
Ehat_ts = Pedgley_data[np.where(Pedgley_data[:,1] == np.max(Pedgley_data[:,1]))[0][0],0]
p_opt = 1,1,1,1,1
delta_ts_opt,_ = spopt.curve_fit(delta_ts, Pedgley_data[:,0], Pedgley_data[:,1], p_opt)
delta_ts_poco_Peg = delta_ts(Ep,*delta_ts_opt)

p_opt = 0,0,1,300,1
delta_e_opt,_ = spopt.curve_fit(delta_e, Farhang_elastic_data[:,0], Farhang_elastic_data[:,1], p_opt, bounds=(0,[1,1,2,500,1]))
delta_e_poco = delta_e(Ep,*delta_e_opt)

p_opt = 0,0,1,300,1,150,1,0.5
delta_re_opt,_ = spopt.curve_fit(delta_re, Pedgley_elastic_data[:,0], Pedgley_elastic_data[:,1], p_opt, bounds=(0,[1,1,2,500,1,500,2,1]))
delta_re_poco = delta_re(Ep,*delta_re_opt)

fig, ax = plt.subplots()
ax.plot(Ep, delta_ts_poco,'-.',c='C0')
ax.plot(Ep, delta_e_poco,'--',c='C0')
ax.plot(Ep, delta_e_poco + delta_ts_poco,c='C0')
ax.plot(Ep, delta_ts_poco_Peg,'-.',c='C2')
ax.plot(Ep, delta_re_poco,'--',c='C2')
ax.plot(Ep, delta_re_poco + delta_ts_poco_Peg,c='C2')
ax.scatter(Farhang_data[:,0],Farhang_data[:,1],c='C0')
ax.scatter(Farhang_elastic_data[:,0],Farhang_elastic_data[:,1],marker='^',c='C0')
ax.scatter(Pedgley_data[:,0],Pedgley_data[:,1],c='C2')
ax.scatter(Pedgley_elastic_data[:,0],Pedgley_elastic_data[:,1],marker='^',c='C2')
ax.grid(True)
ax.set_axisbelow(True)
plt.xlabel('$E_p \ [eV]$')
plt.ylabel('$\delta$')
plt.xlim((0,2000))
plt.ylim((0,0.9))
plt.legend(('$\delta_{ts}$ - Farhang 1993','$\delta_{e}+\delta_{r}$ - Farhang 1993','$\delta$ - Farhang 1993', \
    '$\delta_{ts}$ - Pedgley 1992','$\delta_{e}+\delta_{r}$ - Pedgley 1992','$\delta$ - Pedgley 1992'))
plt.show()

Ee = np.linspace(0,25,100000)

Patino_spectrum_data = np.loadtxt('../data/Graphite/Patino_2015.txt', dtype=(str), delimiter=',').astype(float)

ddeltadE = Ee/(Ee + Phi)**4 # Chung-Everhart fit

fwhm_ind = np.where(np.abs(ddeltadE - np.max(ddeltadE)/2) < 1e-7)[0]
fwhm = Ee[fwhm_ind[1]] - Ee[fwhm_ind[0]]
print(f'FWHM for graphite is {fwhm}')
peak_ind = np.where(ddeltadE == np.max(ddeltadE))[0][0]
print(f'Epeak for graphite is {Ee[peak_ind]}')

plt.plot(Ee,ddeltadE/np.max(Ee/(Ee + Phi)**4))
plt.scatter(Patino_spectrum_data[:,0],Patino_spectrum_data[:,1]/np.max(Patino_spectrum_data[:,1]))
plt.hlines(0.5, xmin=Ee[fwhm_ind[0]], xmax=Ee[fwhm_ind[1]], \
    linestyle='--', color="r", label="FWHM = %0.2f"%fwhm)
plt.axvline(Ee[peak_ind], linestyle='--', color="k", label="$E_{peak}$ = %0.2f"%Ee[peak_ind])
plt.grid()
plt.xlabel('$E_e \ [eV]$')
plt.ylabel('$d\delta / dE_e \ [eV^{-1}]$')
plt.xlim((0,25))
plt.show()