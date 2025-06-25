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
    E = delta_e(Ep,Ehat_e,P1e_inf,P1,W,p)
    R = delta_r(Ep,Er,P1r_inf,r)
    return E+R

Ep = np.linspace(0,5000,100000)
mu = 1

# Tungsten
Phi = 4.55 # Work Function [eV]

Walker_data = np.loadtxt("../data/Tungsten/Walker_2008.txt", dtype=(str), delimiter=',').astype(float)
ElGomati_elastic_data = np.loadtxt("../data/Tungsten/ElGomati_2008.txt", dtype=(str), delimiter=',').astype(float)

deltahat_ts = np.max(Walker_data[:,1])
Ehat_ts = Walker_data[np.where(Walker_data[:,1] == np.max(Walker_data[:,1]))[0][0],0]

p_opt = 1,1,1,1,1
delta_ts_opt,_ = spopt.curve_fit(delta_ts, Walker_data[:,0], Walker_data[:,1], p_opt)
delta_ts_W = delta_ts(Ep,*delta_ts_opt)

p_opt = 0,0,1,300,1,150,1,2
delta_re_opt,_ = spopt.curve_fit(delta_re, ElGomati_elastic_data[:,0], ElGomati_elastic_data[:,1], p_opt, bounds=(0,[1,1,2,500,1,500,2,5]))
# delta_re_W = delta_re(Ep,*delta_re_opt)
delta_re_W = delta_re(Ep,700,0.2,0.06,4900,2,600,0.31,1.2)

fig, ax = plt.subplots()
ax.plot(Ep, delta_ts_W,'-.')
ax.plot(Ep,delta_re_W,'--')
ax.plot(Ep,delta_ts_W+delta_re_W)
ax.scatter(Walker_data[:,0],Walker_data[:,1])
ax.scatter(ElGomati_elastic_data[:,0],ElGomati_elastic_data[:,1])
plt.legend(('$\delta_{ts}$ - Walker 2008', '$\delta_{e} + \delta_{r}$ - El Gomati 2008', '$\delta$'))
ax.grid(True)
ax.set_axisbelow(True)
plt.xlabel('$E_p \ [eV]$')
plt.ylabel('$\delta$')
plt.xlim((0,4000))
plt.ylim((0,1.4))
plt.show()

Ee = np.linspace(0,25,100000)

ddeltadE = Ee/(Ee + Phi)**4 # Chung-Everhart fit

fwhm_ind = np.where(np.abs(ddeltadE - np.max(ddeltadE)/2) < 2e-7)[0]
fwhm = Ee[fwhm_ind[2]] - Ee[fwhm_ind[0]]
print(f'FWHM for tungsten is {fwhm}')
peak_ind = np.where(ddeltadE == np.max(ddeltadE))[0][0]
print(f'Epeak for tungsten is {Ee[peak_ind]}')

plt.plot(Ee,ddeltadE/np.max(Ee/(Ee + Phi)**4))
plt.hlines(0.5, xmin=Ee[fwhm_ind[0]], xmax=Ee[fwhm_ind[2]], \
    linestyle='--', color="r", label="FWHM = %0.2f"%fwhm)
plt.axvline(Ee[peak_ind], linestyle='--', color="k", label="$E_{peak}$ = %0.2f"%Ee[peak_ind])
plt.grid()
plt.xlabel('$E_e \ [eV]$')
plt.ylabel('$d\delta / dE_e \ [eV^{-1}]$')
plt.xlim((0,25))
plt.show()