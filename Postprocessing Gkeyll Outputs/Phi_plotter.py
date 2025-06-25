import numpy as np
import scipy.constants as const
import scipy.integrate as spint
import matplotlib.pyplot as plt
import postgkyl as pg

def load_field(fname):
    data = pg.data.GData(fname)
    dg = pg.GInterpModal(data)
    dg.interpolate(overwrite=True)
    grid = data.get_grid()
    return grid[0][:-1]/lambda_D, data.get_values().reshape(3072,)

phi_b = np.array((0,2,4,6,8,10)) # Bias potentials [kV]
e = const.elementary_charge
eps0 = const.epsilon_0
mi = const.proton_mass
me = const.electron_mass
T0 = 2000*e
n0 = 1.1e23
lambda_D = np.sqrt(eps0*T0/(n0*e**2))
field = {}
plt.rcParams.update({'font.size':18})

# Load No Emission data
for i in phi_b:
    frame = 100
    if i == 6:
        frame = 99
    fname = '/home/lucio/../../mnt/d/NoEmission/NoEmission/' + f'{i}' + 'kV/vp_biasedSheath_NoEmission_' \
        + f'{i}' + 'kV_1x1v_p2-field_' + f'{frame}' + '.gkyl'
    x, field[f'field_{i}kV'] = load_field(fname)


# Load Graphite data
for i in phi_b:
    frame = 100
    if i == 8 or 10:
        frame = 99
    fname = '/home/lucio/../../mnt/d/graphite_H_IIEE/graphite_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-field_' + f'{frame}' + '.gkyl'
    _, field[f'field_g_{i}kV'] = load_field(fname)

# Load Tungsten data
for i in phi_b:
    frame = 100
    fname = '/home/lucio/../../mnt/d/tungsten_H_IIEE/tungsten_H_IIEE/' + f'{i}' + 'kV/vp_biasedSheath_IIEE_' \
        + f'{i}' + 'kV_1x1v_p2-field_' + f'{frame}' + '.gkyl'
    _, field[f'field_w_{i}kV'] = load_field(fname)

phi0 = 1e-3*field['field_0kV']
phi2 = 1e-3*field['field_2kV']
phi4 = 1e-3*field['field_4kV']
phi6 = 1e-3*field['field_6kV']
phi8 = 1e-3*field['field_8kV']
phi10 = 1e-3*field['field_10kV']

phi0_w = 1e-3*field['field_w_0kV']
phi2_w = 1e-3*field['field_w_2kV']
phi4_w = 1e-3*field['field_w_4kV']
phi6_w = 1e-3*field['field_w_6kV']
phi8_w = 1e-3*field['field_w_8kV']
phi10_w = 1e-3*field['field_w_10kV']

phi0_g = 1e-3*field['field_g_0kV']
phi2_g = 1e-3*field['field_g_2kV']
phi4_g = 1e-3*field['field_g_4kV']
phi6_g = 1e-3*field['field_g_6kV']
phi8_g = 1e-3*field['field_g_8kV']
phi10_g = 1e-3*field['field_g_10kV']

plt.plot(x,phi0,'--',c='C0')
plt.scatter(x[phi0==np.max(phi0)],np.max(phi0),c='C0',marker='s',facecolors='none')
plt.plot(x,phi2,'--',c='C1')
plt.scatter(x[phi2==np.max(phi2)],np.max(phi2),c='C1',marker='s',facecolors='none')
plt.plot(x,phi4,'--',c='C2')
plt.scatter(x[phi4==np.max(phi4)],np.max(phi4),c='C2',marker='s',facecolors='none')
plt.plot(x,phi6,'--',c='C3')
plt.scatter(x[phi6==np.max(phi6)],np.max(phi6),c='C3',marker='s',facecolors='none')
plt.plot(x,phi8,'--',c='C4')
plt.scatter(x[phi8==np.max(phi8)],np.max(phi8),c='C4',marker='s',facecolors='none')
plt.plot(x,phi10,'--',c='C5')
plt.scatter(x[phi10==np.max(phi10)],np.max(phi10),c='C5',marker='s',facecolors='none')
plt.plot(x,phi0_w,c='C0',label='0 kV')
plt.scatter(x[phi0_w==np.max(phi0_w)],np.max(phi0_w),c='C0')
plt.plot(x,phi2_w,c='C1',label='2 kV')
plt.scatter(x[phi2_w==np.max(phi2_w)],np.max(phi2_w),c='C1')
plt.plot(x,phi4_w,c='C2',label='4 kV')
plt.scatter(x[phi4_w==np.max(phi4_w)],np.max(phi4_w),c='C2')
plt.plot(x,phi6_w,c='C3',label='6 kV')
plt.scatter(x[phi6_w==np.max(phi6_w)],np.max(phi6_w),c='C3')
plt.plot(x,phi8_w,c='C4',label='8 kV')
plt.scatter(x[phi8_w==np.max(phi8_w)],np.max(phi8_w),c='C4')
plt.plot(x,phi10_w,c='C5',label='10 kV')
plt.scatter(x[phi10_w==np.max(phi10_w)],np.max(phi10_w),c='C5')
plt.plot(x,phi0_g,'-.',c='C0')
plt.scatter(x[phi0_g==np.max(phi0_g)],np.max(phi0_g),c='C0',marker='^',facecolors='none')
plt.plot(x,phi2_g,'-.',c='C1')
plt.scatter(x[phi2_g==np.max(phi2_g)],np.max(phi2_g),c='C1',marker='^',facecolors='none')
plt.plot(x,phi4_g,'-.',c='C2')
plt.scatter(x[phi4_g==np.max(phi4_g)],np.max(phi4_g),c='C2',marker='^',facecolors='none')
plt.plot(x,phi6_g,'-.',c='C3')
plt.scatter(x[phi6_g==np.max(phi6_g)],np.max(phi6_g),c='C3',marker='^',facecolors='none')
plt.plot(x,phi8_g,'-.',c='C4')
plt.scatter(x[phi8_g==np.max(phi8_g)],np.max(phi8_g),c='C4',marker='^',facecolors='none')
plt.plot(x,phi10_g,'-.',c='C5')
plt.scatter(x[phi10_g==np.max(phi10_g)],np.max(phi10_g),c='C5',marker='^',facecolors='none')
plt.xlim([np.min(x),np.max(x)])
plt.ylim([0,14])
plt.xlabel('$x / \lambda_D$')
plt.ylabel('$\phi_p \ [kV]$')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title('Plasma vs Bias Potential')
plt.grid()
plt.show()

fig,ax = plt.subplots(1,2)
ax[1].plot(x,phi8,'--',c='C1',lw=2)
ax[1].plot(x,phi8_w,c='C0',lw=2)
ax[1].plot(x,phi8_g,'-.',c='C2',lw=2)
ax[1].set_xlim([240,np.max(x)])
ax[1].axvline(x[318],ls='--',c='C0',lw=2)
ax[1].set_ylim([0,10])
ax[0].plot(x,phi8,'--',c='C1',lw=2)
ax[0].plot(x,phi8_w,c='C0',lw=2)
ax[0].plot(x,phi8_g,'-.',c='C2',lw=2)
ax[0].axvline(x[247],ls='--',c='C0',lw=2)
ax[0].set_xlim([min(x),-240])
ax[0].set_ylim([8,8.7])
ax[0].set_xlabel('$x / \lambda_D$')
ax[1].set_xlabel('$x / \lambda_D$')
ax[0].set_ylabel('$\phi_p \ [kV]$')
ax[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax[1].legend(('No Emission','Tungsten','Graphite'))
# plt.title('Plasma vs Bias Potential')
ax[0].grid()
ax[1].grid()
plt.show()