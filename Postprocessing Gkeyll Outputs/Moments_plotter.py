import numpy as np
import scipy.constants as const
import scipy.integrate as spint
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import postgkyl as pg


def load_mom(fname):
    data = pg.data.GData(fname)
    dg = pg.GInterpModal(data,num_interp=i_pts)
    dg.interpolate(overwrite=True)
    grid = data.get_grid()
    return grid[0][:-1]/lambda_D, data.get_values().reshape(1024*i_pts,)

def load_distFun(fname):
    data = pg.data.GData(fname)
    dg = pg.GInterpModal(data)
    dg.interpolate(overwrite=True)
    grid = data.get_grid()
    return grid[0][:-1]/lambda_D, grid[1][:-1], data.get_values()[0,:,0]

i_pts = 6
phi_b = np.array((0,2,4,6,8,10)) # Bias potentials [kV]
e = const.elementary_charge
eps0 = const.epsilon_0
mi = const.proton_mass
me = const.electron_mass
T0 = 2000*e
n0 = 1.1e23
lambda_D = np.sqrt(eps0*T0/(n0*e**2))
v_the = np.sqrt(T0/me)
v_thi = np.sqrt(T0/mi)
ion_M0 = {}
elc_M0 = {}
ion_M1 = {}
elc_M1 = {}
ion_M2 = {}
elc_M2 = {}
ion_M3 = {}
elc_M3 = {}
Ti = {}
Te = {}
qi = {}
qe = {}
plt.rcParams.update({'font.size':18})

# Load No Emission data
for i in phi_b:
    frame = 100
    if i == 6:
        frame = 99
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M0_' + f'{frame}' + '.gkyl'
    x, ion_M0[f'{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M0_' + f'{frame}' + '.gkyl'
    _, elc_M0[f'{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    _, ion_M1[f'{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M1i_' + f'{frame}' + '.gkyl'
    _, elc_M1[f'{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M2_' + f'{frame}' + '.gkyl'
    _, ion_M2[f'{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M2_' + f'{frame}' + '.gkyl'
    _, elc_M2[f'{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M3i_' + f'{frame}' + '.gkyl'
    _, ion_M3[f'{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M3i_' + f'{frame}' + '.gkyl'
    _, elc_M3[f'{i}kV'] = load_mom(fname)
    Ti[f'{i}kV'] = mi*(ion_M2[f'{i}kV']/ion_M0[f'{i}kV'] - (ion_M1[f'{i}kV']/ion_M0[f'{i}kV'])**2)
    Te[f'{i}kV'] = me*(elc_M2[f'{i}kV']/elc_M0[f'{i}kV'] - (elc_M1[f'{i}kV']/elc_M0[f'{i}kV'])**2)
    qi[f'{i}kV'] = mi*(0.5*ion_M3[f'{i}kV'] - ion_M1[f'{i}kV']*(1.5*ion_M2[f'{i}kV']/ion_M0[f'{i}kV'] \
        - (ion_M1[f'{i}kV']/ion_M0[f'{i}kV'])**2))
    qe[f'{i}kV'] = me*(0.5*elc_M3[f'{i}kV'] - elc_M1[f'{i}kV']*(1.5*elc_M2[f'{i}kV']/elc_M0[f'{i}kV'] \
        - (elc_M1[f'{i}kV']/elc_M0[f'{i}kV'])**2))

# Load Graphite data
for i in phi_b:
    frame = 100
    if i == 8 or 10:
        frame = 99
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M0_' + f'{frame}' + '.gkyl'
    _, ion_M0[f'g_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M0_' + f'{frame}' + '.gkyl'
    _, elc_M0[f'g_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    _, ion_M1[f'g_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M1i_' + f'{frame}' + '.gkyl'
    _, elc_M1[f'g_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M2_' + f'{frame}' + '.gkyl'
    _, ion_M2[f'g_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M2_' + f'{frame}' + '.gkyl'
    _, elc_M2[f'g_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M3i_' + f'{frame}' + '.gkyl'
    _, ion_M3[f'g_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M3i_' + f'{frame}' + '.gkyl'
    _, elc_M3[f'g_{i}kV'] = load_mom(fname)
    Ti[f'g_{i}kV'] = mi*(ion_M2[f'g_{i}kV']/ion_M0[f'g_{i}kV'] - (ion_M1[f'g_{i}kV']/ion_M0[f'g_{i}kV'])**2)
    Te[f'g_{i}kV'] = me*(elc_M2[f'g_{i}kV']/elc_M0[f'g_{i}kV'] - (elc_M1[f'g_{i}kV']/elc_M0[f'g_{i}kV'])**2)
    qi[f'g_{i}kV'] = mi*(0.5*ion_M3[f'g_{i}kV'] - ion_M1[f'g_{i}kV']*(1.5*ion_M2[f'g_{i}kV']/ion_M0[f'g_{i}kV'] \
        - (ion_M1[f'g_{i}kV']/ion_M0[f'g_{i}kV'])**2))
    qe[f'g_{i}kV'] = me*(0.5*elc_M3[f'g_{i}kV'] - elc_M1[f'g_{i}kV']*(1.5*elc_M2[f'g_{i}kV']/elc_M0[f'g_{i}kV'] \
        - (elc_M1[f'g_{i}kV']/elc_M0[f'g_{i}kV'])**2))

# Load Tungsten data
for i in phi_b:
    frame = 100
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M0_' + f'{frame}' + '.gkyl'
    _, ion_M0[f'w_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M0_' + f'{frame}' + '.gkyl'
    _, elc_M0[f'w_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M1i_' + f'{frame}' + '.gkyl'
    _, ion_M1[f'w_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M1i_' + f'{frame}' + '.gkyl'
    _, elc_M1[f'w_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M2_' + f'{frame}' + '.gkyl'
    _, ion_M2[f'w_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M2_' + f'{frame}' + '.gkyl'
    _, elc_M2[f'w_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-ion_M3i_' + f'{frame}' + '.gkyl'
    _, ion_M3[f'w_{i}kV'] = load_mom(fname)
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-elc_M3i_' + f'{frame}' + '.gkyl'
    _, elc_M3[f'w_{i}kV'] = load_mom(fname)
    Ti[f'w_{i}kV'] = mi*(ion_M2[f'w_{i}kV']/ion_M0[f'w_{i}kV'] - (ion_M1[f'w_{i}kV']/ion_M0[f'w_{i}kV'])**2)
    Te[f'w_{i}kV'] = me*(elc_M2[f'w_{i}kV']/elc_M0[f'w_{i}kV'] - (elc_M1[f'w_{i}kV']/elc_M0[f'w_{i}kV'])**2)
    qi[f'w_{i}kV'] = mi*(0.5*ion_M3[f'w_{i}kV'] - ion_M1[f'w_{i}kV']*(1.5*ion_M2[f'w_{i}kV']/ion_M0[f'w_{i}kV'] \
        - (ion_M1[f'w_{i}kV']/ion_M0[f'w_{i}kV'])**2))
    qe[f'w_{i}kV'] = me*(0.5*elc_M3[f'w_{i}kV'] - elc_M1[f'w_{i}kV']*(1.5*elc_M2[f'w_{i}kV']/elc_M0[f'w_{i}kV'] \
        - (elc_M1[f'w_{i}kV']/elc_M0[f'w_{i}kV'])**2))


### Density ###

fig,ax = plt.subplots(3,2,figsize=(12,10))
k=0
for i in phi_b[[0,2,-2]]:
    ax[k,0].plot(x,ion_M0[f'{i}kV']/n0,'--',c='b')
    ax[k,0].plot(x,ion_M0[f'w_{i}kV']/n0,c='b')
    ax[k,0].plot(x,ion_M0[f'g_{i}kV']/n0,'-.',c='b')
    ax[k,0].plot(x,elc_M0[f'{i}kV']/n0,'--',c='r')
    ax[k,0].plot(x,elc_M0[f'w_{i}kV']/n0,c='r')
    ax[k,0].plot(x,elc_M0[f'g_{i}kV']/n0,'-.',c='r')
    ax[k,1].plot(x,ion_M0[f'{i}kV']/n0,'--',c='b')
    ax[k,1].plot(x,ion_M0[f'w_{i}kV']/n0,c='b')
    ax[k,1].plot(x,ion_M0[f'g_{i}kV']/n0,'-.',c='b')
    ax[k,1].plot(x,elc_M0[f'{i}kV']/n0,'--',c='r')
    ax[k,1].plot(x,elc_M0[f'w_{i}kV']/n0,c='r')
    ax[k,1].plot(x,elc_M0[f'g_{i}kV']/n0,'-.',c='r')
    k += 1
for a in ax[:, 0]:
    a.set_xlim((-256, -240))
    a.set_ylim((0, 0.65))
    a.set_ylabel('$n/n_0$')
    a.grid()
for a in ax[:, 1]:
    a.set_xlim((240, 256))
    a.set_ylim((0, 0.65))
    a.set_ylabel('$n/n_0$')
    a.grid()
for a in ax[2, :]:
    a.set_xlabel('$x/\lambda_D$')
plt.tight_layout()
plt.show()

### Temperature ###

# Electrons
fig,ax = plt.subplots(1,2,figsize=(16,6))
k = 0
for i in phi_b[[0,2,-1]]:
    ax[0].plot(x,Te[f'{i}kV']/T0,'--',c=f'C{k}')
    ax[0].plot(x,Te[f'w_{i}kV']/T0,c=f'C{k}',label=f'{i} kV')
    ax[0].plot(x,Te[f'g_{i}kV']/T0,'-.',c=f'C{k}')
    ax[1].plot(x,Te[f'{i}kV']/T0,'--',c=f'C{k}')
    ax[1].plot(x,Te[f'w_{i}kV']/T0,c=f'C{k}',label=f'{i} kV')
    ax[1].plot(x,Te[f'g_{i}kV']/T0,'-.',c=f'C{k}')
    k += 1
ax[0].set_xlabel('$x/\lambda_D$')
ax[1].set_xlabel('$x/\lambda_D$')
ax[0].set_ylabel('$T_e/T_0$')
ax[1].set_ylabel('$T_e/T_0$')
ax[0].set_xlim((-256,-240))
ax[0].set_ylim((0,0.4))
ax[1].set_xlim((240,256))
ax[1].set_ylim((0,2))
ax[0].legend(loc='upper left')
ax[1].legend(loc='upper left')
ax[0].grid()
ax[1].grid()
plt.show()

# Ions
fig,ax = plt.subplots(1,2,figsize=(16,6))
k = 0
for i in phi_b[[0,2,-1]]:
    ax[0].plot(x,Ti[f'{i}kV']/T0,'--',c=f'C{k}')
    ax[0].plot(x,Ti[f'w_{i}kV']/T0,c=f'C{k}',label=f'{i} kV')
    ax[0].plot(x,Ti[f'g_{i}kV']/T0,'-.',c=f'C{k}')
    ax[1].plot(x,Ti[f'{i}kV']/T0,'--',c=f'C{k}')
    ax[1].plot(x,Ti[f'w_{i}kV']/T0,c=f'C{k}',label=f'{i} kV')
    ax[1].plot(x,Ti[f'g_{i}kV']/T0,'-.',c=f'C{k}')
    k += 1
ax[0].set_xlabel('$x/\lambda_D$')
ax[1].set_xlabel('$x/\lambda_D$')
ax[0].set_ylabel('$T_i/T_0$')
ax[1].set_ylabel('$T_i/T_0$')
ax[0].set_xlim((-256,-240))
ax[0].set_ylim((0.18,0.3))
ax[1].set_xlim((240,256))
ax[1].set_ylim((0.02,0.3))
ax[0].legend()
ax[1].legend()
ax[0].grid()
ax[1].grid()
plt.show()

### Particle Flux ###

fig,ax = plt.subplots(3,2,figsize=(12,10))
k=0
for i in phi_b[[0,2,-2]]:
    ax[k,0].plot(x,ion_M1[f'{i}kV']/n0/v_the,'--',c='b')
    ax[k,0].plot(x,ion_M1[f'w_{i}kV']/n0/v_the,c='b')
    ax[k,0].plot(x,ion_M1[f'g_{i}kV']/n0/v_the,'-.',c='b')
    ax[k,0].plot(x,elc_M1[f'{i}kV']/n0/v_the,'--',c='r')
    ax[k,0].plot(x,elc_M1[f'w_{i}kV']/n0/v_the,c='r')
    ax[k,0].plot(x,elc_M1[f'g_{i}kV']/n0/v_the,'-.',c='r')
    ax[k,1].plot(x,ion_M1[f'{i}kV']/n0/v_the,'--',c='b')
    ax[k,1].plot(x,ion_M1[f'w_{i}kV']/n0/v_the,c='b')
    ax[k,1].plot(x,ion_M1[f'g_{i}kV']/n0/v_the,'-.',c='b')
    ax[k,1].plot(x,elc_M1[f'{i}kV']/n0/v_the,'--',c='r')
    ax[k,1].plot(x,elc_M1[f'w_{i}kV']/n0/v_the,c='r')
    ax[k,1].plot(x,elc_M1[f'g_{i}kV']/n0/v_the,'-.',c='r')
    k += 1
for a in ax[:, 0]:
    a.set_xlim((-256, -240))
    a.set_ylabel(r'$\Gamma_{\alpha}/n_0 v_{te}$')
    a.grid()
for a in ax[:, 1]:
    a.set_xlim((240, 256))
    a.set_ylabel(r'$\Gamma_{\alpha}/n_0 v_{te}$')
    a.grid()
for a in ax[2, :]:
    a.set_xlabel('$x/\lambda_D$')
plt.tight_layout()
plt.show()

### Heat Flux ###

# Electrons
fig,ax = plt.subplots(3,2,figsize=(12,10))
k=0
for i in phi_b[[0,2,-2]]:
    ax[k,0].plot(x,qe[f'{i}kV'],'--',c='r')
    ax[k,0].plot(x,qe[f'w_{i}kV'],c='r')
    ax[k,0].plot(x,qe[f'g_{i}kV'],'-.',c='r')
    ax[k,1].plot(x,qe[f'{i}kV'],'--',c='r')
    ax[k,1].plot(x,qe[f'w_{i}kV'],c='r')
    ax[k,1].plot(x,qe[f'g_{i}kV'],'-.',c='r')
    k += 1
for a in ax[:, 0]:
    a.set_xlim((-256, -240))
    a.set_ylabel('$q_e \ [W/m^2]$')
    a.grid()
for a in ax[:, 1]:
    a.set_xlim((240, 256))
    a.set_ylabel('$q_e \ [W/m^2]$')
    a.grid()
for a in ax[2, :]:
    a.set_xlabel('$x/\lambda_D$')
plt.tight_layout()
plt.show()

# fig,ax = plt.subplots(1,2,figsize=(16,6))
# ax[0].plot(x,qe['0kV'],'--',c='C0',lw=2)
# ax[0].plot(x,qe['w_0kV'],c='C1',lw=2)
# ax[0].plot(x,qe['g_0kV'],'-.',c='C2',lw=2)
# ax[1].plot(x,qe['0kV'],'--',c='C0',lw=2)
# ax[1].plot(x,qe['w_0kV'],c='C1',lw=2)
# ax[1].plot(x,qe['g_0kV'],'-.',c='C2',lw=2)
# ax[0].set_xlim((-256, -240))
# ax[0].set_ylabel('$q_e \ [W/m^2]$')
# ax[0].grid()
# ax[1].set_xlim((240, 256))
# ax[1].set_ylabel('$q_e \ [W/m^2]$')
# ax[1].grid()
# ax[0].set_xlabel('$x/\lambda_D$')
# ax[1].set_xlabel('$x/\lambda_D$')
# ax[1].legend(('No Emission','Tungsten','Graphite'))
# ax[0].text(0.01, 0.95, '0 kV, anode', transform=ax[0].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# ax[1].text(0.01, 0.95, '0 kV, cathode', transform=ax[1].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# plt.tight_layout()
# plt.show()

# fig,ax = plt.subplots(1,2,figsize=(16,6))
# ax[0].plot(x,qe['4kV'],'--',c='C0',lw=2)
# ax[0].plot(x,qe['w_4kV'],c='C1',lw=2)
# ax[0].plot(x,qe['g_4kV'],'-.',c='C2',lw=2)
# ax[1].plot(x,qe['4kV'],'--',c='C0',lw=2)
# ax[1].plot(x,qe['w_4kV'],c='C1',lw=2)
# ax[1].plot(x,qe['g_4kV'],'-.',c='C2',lw=2)
# ax[0].set_xlim((-256, -240))
# ax[0].set_ylabel('$q_e \ [W/m^2]$')
# ax[0].grid()
# ax[1].set_xlim((240, 256))
# ax[1].set_ylabel('$q_e \ [W/m^2]$')
# ax[1].grid()
# ax[0].set_xlabel('$x/\lambda_D$')
# ax[1].set_xlabel('$x/\lambda_D$')
# ax[1].legend(('No Emission','Tungsten','Graphite'))
# ax[0].text(0.01, 0.95, '4 kV, anode', transform=ax[0].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# ax[1].text(0.01, 0.95, '4 kV, cathode', transform=ax[1].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# plt.tight_layout()
# plt.show()

# fig,ax = plt.subplots(1,2,figsize=(16,6))
# ax[0].plot(x,qe['8kV'],'--',c='C0',lw=2)
# ax[0].plot(x,qe['w_8kV'],c='C1',lw=2)
# ax[0].plot(x,qe['g_8kV'],'-.',c='C2',lw=2)
# ax[1].plot(x,qe['8kV'],'--',c='C0',lw=2)
# ax[1].plot(x,qe['w_8kV'],c='C1',lw=2)
# ax[1].plot(x,qe['g_8kV'],'-.',c='C2',lw=2)
# ax[0].set_xlim((-256, -240))
# ax[0].set_ylabel('$q_e \ [W/m^2]$')
# ax[0].grid()
# ax[1].set_xlim((240, 256))
# ax[1].set_ylabel('$q_e \ [W/m^2]$')
# ax[1].grid()
# ax[0].set_xlabel('$x/\lambda_D$')
# ax[1].set_xlabel('$x/\lambda_D$')
# ax[1].legend(('No Emission','Tungsten','Graphite'))
# ax[0].text(0.01, 0.95, '8 kV, anode', transform=ax[0].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# ax[1].text(0.01, 0.95, '8 kV, cathode', transform=ax[1].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# plt.tight_layout()
# plt.show()

# Ions
fig,ax = plt.subplots(3,2,figsize=(12,10))
k=0
for i in phi_b[[0,2,-2]]:
    ax[k,0].plot(x,qi[f'{i}kV'],'--',c='b')
    ax[k,0].plot(x,qi[f'w_{i}kV'],c='b')
    ax[k,0].plot(x,qi[f'g_{i}kV'],'-.',c='b')
    ax[k,1].plot(x,qi[f'{i}kV'],'--',c='b')
    ax[k,1].plot(x,qi[f'w_{i}kV'],c='b')
    ax[k,1].plot(x,qi[f'g_{i}kV'],'-.',c='b')
    k += 1
for a in ax[:, 0]:
    a.set_xlim((-256, -240))
    a.set_ylabel('$q_i \ [W/m^2]$')
    a.grid()
for a in ax[:, 1]:
    a.set_xlim((240, 256))
    a.set_ylabel('$q_i \ [W/m^2]$')
    a.grid()
for a in ax[2, :]:
    a.set_xlabel('$x/\lambda_D$')
plt.tight_layout()
plt.show()

# fig,ax = plt.subplots(1,2,figsize=(16,6))
# ax[0].plot(x,qi['0kV'],'--',c='C0',lw=2)
# ax[0].plot(x,qi['w_0kV'],c='C1',lw=2)
# ax[0].plot(x,qi['g_0kV'],'-.',c='C2',lw=2)
# ax[1].plot(x,qi['0kV'],'--',c='C0',lw=2)
# ax[1].plot(x,qi['w_0kV'],c='C1',lw=2)
# ax[1].plot(x,qi['g_0kV'],'-.',c='C2',lw=2)
# ax[0].set_xlim((-256, -240))
# ax[0].set_ylabel('$q_i \ [W/m^2]$')
# ax[0].grid()
# ax[1].set_xlim((240, 256))
# ax[1].set_ylabel('$q_i \ [W/m^2]$')
# ax[1].grid()
# ax[0].set_xlabel('$x/\lambda_D$')
# ax[1].set_xlabel('$x/\lambda_D$')
# ax[1].legend(('No Emission','Tungsten','Graphite'))
# ax[0].text(0.01, 0.95, '0 kV, anode', transform=ax[0].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# ax[1].text(0.01, 0.95, '0 kV, cathode', transform=ax[1].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# plt.tight_layout()
# plt.show()

# fig,ax = plt.subplots(1,2,figsize=(16,6))
# ax[0].plot(x,qi['4kV'],'--',c='C0',lw=2)
# ax[0].plot(x,qi['w_4kV'],c='C1',lw=2)
# ax[0].plot(x,qi['g_4kV'],'-.',c='C2',lw=2)
# ax[1].plot(x,qi['4kV'],'--',c='C0',lw=2)
# ax[1].plot(x,qi['w_4kV'],c='C1',lw=2)
# ax[1].plot(x,qi['g_4kV'],'-.',c='C2',lw=2)
# ax[0].set_xlim((-256, -240))
# ax[0].set_ylabel('$q_i \ [W/m^2]$')
# ax[0].grid()
# ax[1].set_xlim((240, 256))
# ax[1].set_ylabel('$q_i \ [W/m^2]$')
# ax[1].grid()
# ax[0].set_xlabel('$x/\lambda_D$')
# ax[1].set_xlabel('$x/\lambda_D$')
# ax[1].legend(('No Emission','Tungsten','Graphite'))
# ax[0].text(0.01, 0.95, '4 kV, anode', transform=ax[0].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# ax[1].text(0.01, 0.95, '4 kV, cathode', transform=ax[1].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# plt.tight_layout()
# plt.show()

# fig,ax = plt.subplots(1,2,figsize=(16,6))
# ax[0].plot(x,qi['8kV'],'--',c='C0',lw=2)
# ax[0].plot(x,qi['w_8kV'],c='C1',lw=2)
# ax[0].plot(x,qi['g_8kV'],'-.',c='C2',lw=2)
# ax[1].plot(x,qi['8kV'],'--',c='C0',lw=2)
# ax[1].plot(x,qi['w_8kV'],c='C1',lw=2)
# ax[1].plot(x,qi['g_8kV'],'-.',c='C2',lw=2)
# ax[0].set_xlim((-256, -240))
# ax[0].set_ylabel('$q_i \ [W/m^2]$')
# ax[0].grid()
# ax[1].set_xlim((240, 256))
# ax[1].set_ylabel('$q_i \ [W/m^2]$')
# ax[1].grid()
# ax[0].set_xlabel('$x/\lambda_D$')
# ax[1].set_xlabel('$x/\lambda_D$')
# ax[1].legend(('No Emission','Tungsten','Graphite'))
# ax[0].text(0.01, 0.95, '8 kV, anode', transform=ax[0].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# ax[1].text(0.01, 0.95, '8 kV, cathode', transform=ax[1].transAxes, fontsize=12,
# verticalalignment='top', horizontalalignment='left',bbox=dict(facecolor='lightgray',\
# edgecolor='none', boxstyle='round,pad=0.2'))
# plt.tight_layout()
# plt.show()