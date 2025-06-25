import numpy as np
import math
import scipy.constants as const
import scipy.optimize as spopt
import scipy.integrate as spint
from scipy.special import digamma, roots_legendre

###### Wall Material Constants ######
e0 = const.epsilon_0
e = const.elementary_charge
Phi = 4.62 # Work function, eV, DOI: https://doi.org/10.1098/rspa.1952.0116
W = Phi # Barrier height, eV
I = 78 # Mean excitation energy, eV, DOI: https://doi.org/10.1093/jicru_os25.2.6
Z = 6 # Target charge state
ma = 12.0107 # Atomic weight [amu]
U = 1.66e-27 # Atomic mass unit [kg]
rho = 2260 # Mass density [kg/m^3], 1780 for POCO AXF 5Q graphite
nw = rho/(ma*U) # Atomic number density [m^-3]
mi = const.proton_mass # Impacting ion mass [kg]

###### SRIM-Schou-Bethe Model ######

# Digamma function
def phi(x):
    return digamma(x)

# Modified Bethe Formula (DOI: https://doi.org/10.1016/j.ultramic.2014.11.003)
def SeFun(Ee):
    E = Ee + W
    logArg = (0.5*np.e)**0.5 * E/I
    C = (2*math.pi*(e**4)*Z*nw)/(E*e*(4*math.pi*e0)**2)
    G = 1 - ((np.e/2)**(0.5))*np.log(1 + (E/I)**2)*(I/E) \
    + (1/3)*np.log(Z/2)*np.exp(-(3/Z**(0.5))*(1 - 2/Z**(0.5) + np.log(E/I))**2)*(E/I)
    return C*np.log(logArg+G)

# Load SRIM stopping power [MeV/mm]
fname = "../data/Graphite/Hydrogen in Graphite carbon (ICRU-906).txt"
data = np.loadtxt(fname, dtype=(str))
data = data[:, (0,2,3)]
data = data.astype(float)

Si = 1.602e-10*data[:,1] # Electronic stopping power [J/m]
Si_n = 1.602e-10*data[:,2] # Nuclear stopping power [J/m]
Ei = data[:, 0] # Ion impact energy
# Convert to keV
eVindex = np.where(Ei == 999.999)
Ei[:eVindex[0][0]+1] = 1e-3*Ei[:eVindex[0][0]+1]
MeVindex = np.where(Ei == 1)
Ei[MeVindex[0][0]:] = 1e3*Ei[MeVindex[0][0]:]

# Calculate wall term, find peak energy
def wallCalc(Ee):
    return Ee/(4*((Ee+W)**2)*SeFun(Ee))
def gaussian(Es, E0, tau):
    return np.exp(-(np.log(Es/E0))**2/(2*tau**2))
Ee = np.logspace(-4,4,100000)
wallTerm = wallCalc(Ee)
wallTerm_Epeak = Ee[int(np.argmax(wallTerm))]
m = 2 - 0.5*(W/wallTerm_Epeak)
Gm = m/(phi(1) - phi(1-m))

# Binding Energy Correction
data = np.loadtxt('../data/Graphite/graphite_DOS.csv',dtype=(str),delimiter=',')
data = data[:, (0,1)]
data = data.astype(float) # Density of States [eV^-1] vs E-Ef [eV], DOI: https://doi.org/10.1088/0953-8984/21/33/335502
E = data[:,0]
DOS = data[:,1]
dE = E[1] - E[0]
def FermiFun(E):
    return 1/(np.exp(E/Tf) + 1)
Tf = 0.025 # 0.025 < Tf < 0.334 [eV], room temp -> sublimation temp
V = spint.trapezoid(-E*DOS*FermiFun(E),dx=dE)/spint.trapezoid(DOS*FermiFun(E),dx=dE) # Average Valence Electron Binding Energy [eV]

def wallCalc(Ee):
    return (Gm*Ee)/(4*(Ee + W)*(Ee + W + (2-m)*V)*SeFun(Ee))
wallTerm = wallCalc(Ee) # Calculate wall term with binding correction

dgammadE = (Si[Ei==500]+Si_n[Ei==500])*wallTerm # Ion-induced spectrum for ion impact energy = 500 keV

# Integrating wall term, Lorentzian fit w/ Gauss-Legendre quadrature
def lorentzian(Es, E0, alpha, beta, tau): # Vector input
    L = (1 /(1 + (np.log(Es/E0))**2/(2*tau**2)))
    L[Es <= E0] = L[Es <= E0]**alpha
    L[Es > E0] = L[Es > E0]**beta
    return L
def lorentzian2(Es, E0, alpha, beta, tau): # Scalar input
    L = (1 /(1 + ((np.log(Es/E0))**2 / (2*tau**2))))
    if Es <= E0:
        L = L**alpha
    else:
        L = L**beta
    return L
p_lorentz = Ee[wallTerm == np.max(wallTerm)][0], 1.0, 1.0, 1.0    # Initial guesses
lorentz_opt,_ = spopt.curve_fit(lorentzian, Ee, wallTerm/np.max(wallTerm), p_lorentz)  # Fit

# Gaussian quadrature integration
nodes,weights = roots_legendre(64)
def gauss_legendre_integrate(func, a, b):
    integral = 0.0
    for i in range(len(nodes)):
        x = 0.5 * (b - a)*nodes[i] + 0.5*(b + a) # Transform the node to the integration interval
        integral += weights[i] * func(x)
    integral *= 0.5*(b - a)  # Scale the result by the interval length
    return integral

g = gauss_legendre_integrate(lambda t: lorentzian2(1/(1-t), \
    *lorentz_opt) * 1/(1-t)**2, 0, 1)*np.max(wallTerm)

gamma = (Si+Si_n)*g # Ion-induced yield w/ SRIM ion stopping power

###### Lindhard-Bethe-Schou Model ######

def SiFun(E): # Ion Stopping power using Modified Bethe Formula
    logArg = (0.5*np.e)**0.5 * E/I
    C = (4*math.pi*(e**4)*Z*nw)/(E*e*(4*math.pi*e0)**2)
    G = 1 - ((np.e/2)**(0.5))*np.log(1 + (E/I)**2)*(I/E) \
    + (1/3)*np.log(Z/2)*np.exp(-(3/Z**(0.5))*(1 - 2/Z**(0.5) + np.log(E/I))**2)*(E/I)
    return C*np.log(logArg+G)

Ei2_eV = np.logspace(0,8,208)
Ei2 = np.logspace(-3,5,208)
vi = np.sqrt(2*Ei2_eV*e / (mi))
T = vi**2 * const.electron_mass / e
Si_b = SiFun(T) # Bethe ion stopping power [J/m]

# Lindhard stopping power, Iafrate and Ziegler (DOI: https://doi.org/10.1063/1.326750)
me = 9.109e-31 * 1e3  # Electron mass [g]
hbar = 6.626e-34 / (2*np.pi) * 1e7 # Reduced Planck const. [ergs*s]
e = 4.80326e-10 # Elementary charge [esu]
U = U * 1e3 # [g/amu]
rho = rho * 1e-3 # [g/cm^3]
nw = rho/(ma*U) # [cm^-3]

# Calculate Fermi velocity
vf = (hbar/me) * (3*nw*np.pi**2)**(1/3) # Fermi velocity [cm/s]
X2 = (e**2) / (hbar * np.pi * vf)

def f1_calc(z, u):
    log_term1 = np.log(np.abs((z - u + 1) / (z - u - 1) ))
    log_term2 = np.log(np.abs((z + u + 1) / (z + u - 1) ))
    f1 = 0.5 + (1/(8*z)) * (1 - (z-u)**2) * log_term1 \
        + (1/(8*z)) * (1 - (z + u)**2) * log_term2
    return f1

def f2_calc(z,u):          
    if (z + u) < 1:
        return np.pi*u/2
    elif np.abs(z-u) < 1 and (z+u) > 1 and np.abs(z-u) < (z+u):
        return np.pi*(1-(z-u)**2)/(8*z)
    elif np.abs(z-u) > 1:
        return 0

def inner_integrand(t, u):
    z = t/(1-t) # Change of variables to avoid infinite limit
    dzdt = 1/(1-t)**2 # Jacobian 
    f1 = f1_calc(z, u)
    f2 = f2_calc(z, u)
    numer = z**3 * f2
    denom = (z**2 + X2*f1)**2 + (X2*f2)**2
    return numer/denom * dzdt

def outer_integrand(u):
    t_min, t_max = 0, 1
    inner_integral_val = gauss_legendre_integrate(lambda t: inner_integrand(t, u), t_min, t_max)
    return u*inner_integral_val

def S_calc(v):
    S = []
    for vel in v:
        u_min, u_max = 0, vel/vf
        L = gauss_legendre_integrate(outer_integrand, u_min, u_max)
        S_val = (4*np.pi / me) * (e**2 / vel)**2 * nw * (6 / np.pi) * L
        S.append(S_val)
    return np.array(S)

v = np.sqrt(2*1e3*Ei2*1.602e-12 / (mi*1e3)) # ion velocity [cm/s]
Si_L = S_calc(v) * 1e-5 # Lindhard ion stopping power [J/m]

gamma_LBS = (Si_L+Si_b+Si_n)*g # Ion-induced yield