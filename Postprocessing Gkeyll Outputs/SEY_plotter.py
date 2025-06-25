import numpy as np
import scipy.constants as const
import scipy.integrate as spint
import matplotlib.pyplot as plt
import postgkyl as pg

def load_flux(fname):
    data = pg.data.GData(fname)
    dg = pg.GInterpModal(data)
    dg.interpolate(overwrite=True)
    return data.get_values()

def load_distFun(fname):
    data = pg.data.GData(fname)
    dg = pg.GInterpModal(data)
    dg.interpolate(overwrite=True)
    grid = data.get_grid()
    return grid[0][:-1], grid[1][:-1], data.get_values()[0,:,0]

phi_b = np.array((0,2,4,6,8,10)) # Bias potentials [kV]
SEY_g_c = []
SEY_g_a = []
SEY_w_c = []
SEY_w_a = []

# Load Graphite data
for i in phi_b:
    frame = 100
    if i == 8 or 10:
        frame = 99
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    ion_flux = load_flux(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_bc_up_' + f'{frame}' + '.gkyl'
    _, V, elc_bc_cat = load_distFun(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_bc_lo_' + f'{frame}' + '.gkyl'
    _, _, elc_bc_ano = load_distFun(fname)
    dv = V[1] - V[0]
    
    elc_flux_ano, elc_flux_cat = spint.trapezoid(elc_bc_ano[V>=0]*V[V>=0],dx=dv), spint.trapezoid(elc_bc_cat[V<=0]*V[V<=0],dx=dv)
    SEY_g_c.append(np.abs(elc_flux_cat/ion_flux[-1][0]))
    SEY_g_a.append(np.abs(elc_flux_ano/ion_flux[0][0]))

# Load Tungsten data
for i in phi_b:
    frame = 100
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    ion_flux = load_flux(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_bc_up_' + f'{frame}' + '.gkyl'
    _, _, elc_bc_cat = load_distFun(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_bc_lo_' + f'{frame}' + '.gkyl'
    _, _, elc_bc_ano = load_distFun(fname)
    dv = V[1] - V[0]
    elc_flux_ano, elc_flux_cat = spint.trapezoid(elc_bc_ano[V>=0]*V[V>=0],dx=dv), spint.trapezoid(elc_bc_cat[V<=0]*V[V<=0],dx=dv)
    SEY_w_c.append(np.abs(elc_flux_cat/ion_flux[-1][0]))
    SEY_w_a.append(np.abs(elc_flux_ano/ion_flux[0][0]))

afit_g = np.polyfit(phi_b,SEY_g_a,deg=1)
cfit_g = np.polyfit(phi_b,SEY_g_c,deg=1)
afit_w = np.polyfit(phi_b,SEY_w_a,deg=1)
cfit_w = np.polyfit(phi_b,SEY_w_c,deg=1)

plt.scatter(phi_b,SEY_w_a)
plt.scatter(phi_b,SEY_w_c)
plt.scatter(phi_b,SEY_g_a)
plt.scatter(phi_b,SEY_g_c)
plt.plot(phi_b,afit_w[0]*phi_b + afit_w[1],'--')
plt.plot(phi_b,cfit_w[0]*phi_b + cfit_w[1],'--')
plt.plot(phi_b,afit_g[0]*phi_b + afit_g[1],'--')
plt.plot(phi_b,cfit_g[0]*phi_b + cfit_g[1],'--')
plt.legend(('Anode - W', 'Cathode - W','Anode - C','Cathode - C'))
plt.xlabel('$\phi_b \ [kV]$')
plt.ylabel('$\gamma$')
# plt.title('Ion-induced Electron Yield at Anode, Cathode')
plt.grid()
plt.show()