import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import postgkyl as pg

def load_flux(fname):
    data = pg.data.GData(fname)
    dg = pg.GInterpModal(data)
    dg.interpolate(overwrite=True)  
    return data.get_values()

def p_0(me, mi, Te, Ti, g, d):
    return -0.5*np.log((2*np.pi*me/mi)*(1 + Ti/Te)*(((1+g)/(1-d))**2))

def p(pb, p0):
    return -np.log(2*np.exp(-p0)/(1 + np.exp(e*pb/Te)))

e = const.elementary_charge # [C]
a = 3e-3 # Pinch radius [m]
me = 9.109e-31
mi = 1.67e-27
Te = 2000*e
Ti = Te
cs = np.sqrt(2*Te/1.67e-27)
phi_b = np.array((0,2,4,6,8,10)) # Bias potentials [kV]
I = []
I_g = []
I_w = []

# Load No Emission data
for i in phi_b:
    frame = 100
    if i == 6:
        frame = 99
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    ion_flux = load_flux(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M1i_' + f'{frame}' + '.gkyl'
    elc_flux = load_flux(fname)
    I.append(1e-3*(np.pi*a**2)*e*np.mean(ion_flux-elc_flux))

# Load Graphite data
for i in phi_b:
    frame = 100
    if i == 8 or 10:
        frame = 99
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    ion_flux = load_flux(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M1i_' + f'{frame}' + '.gkyl'
    elc_flux = load_flux(fname)
    I_g.append(1e-3*(np.pi*a**2)*e*np.mean(ion_flux-elc_flux))

# Load Tungsten data
for i in phi_b:
    frame = 100
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    ion_flux = load_flux(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' + f'{i}' + 'kV_1x1v_p2-elc_M1i_' + f'{frame}' + '.gkyl'
    elc_flux = load_flux(fname)
    I_w.append(1e-3*(np.pi*a**2)*e*np.mean(ion_flux-elc_flux))

phi_b0 = np.linspace(0,10,100)
gamma_w = 0.05541898835650642*phi_b0
gamma_g = 0.13866926910245958*phi_b0

# Theoretical Model
I0 = 1e-3*0.5*np.pi*(0.003**2)*e*(1.1e23)*cs* \
    (1-np.exp(p_0(me,mi,Te,Ti,0,0) - (p(1e3*phi_b0,p_0(me,mi,Te,Ti,0,0)))))
I0_g = 1e-3*0.5*np.pi*(0.003**2)*e*(1.1e23)*cs* \
    ((1+gamma_g)-np.exp(p_0(me,mi,Te,Ti,gamma_g,0) - (p(1e3*phi_b0,p_0(me,mi,Te,Ti,gamma_g,0)))))
I0_w = 1e-3*0.5*np.pi*(0.003**2)*e*(1.1e23)*cs* \
    ((1+gamma_w)-np.exp(p_0(me,mi,Te,Ti,gamma_w,0) - (p(1e3*phi_b0,p_0(me,mi,Te,Ti,gamma_w,0)))))

plt.plot(phi_b0,I0,'--')
plt.plot(phi_b0,I0_w,'--',label='_nolegend_')
plt.plot(phi_b0,I0_g,'--',label='_nolegend_')
plt.scatter(phi_b,I, zorder=2) 
plt.scatter(phi_b,I_w, zorder=2)
plt.scatter(phi_b,I_g, zorder=2)
plt.errorbar(4,125,35,marker='d',mec='k',c='r',ecolor='k',capsize=5, capthick=1, zorder=2)
plt.errorbar(5,250,50,marker='d',mec='k',c='r',ecolor='k',capsize=5, capthick=1, zorder=2)
plt.scatter(6,200,marker='d',edgecolors='k',c='r', zorder=2)
plt.legend(('Theory','No Emission','Tungsten','Graphite','Experiment'))
plt.xlabel('$\phi_b \ [kV]$')
plt.ylabel('$I \ [kA]$')
# plt.title('Pinch Current vs Bias Potential')
plt.xlim((-0.1,10.1))
plt.grid(zorder=1)
plt.show()