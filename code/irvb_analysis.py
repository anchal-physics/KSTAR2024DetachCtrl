import h5py
import numpy as np
import mat4py
import scipy
import matplotlib
import matplotlib.pyplot as plt
import shapely
import shapely.geometry
import copy

efit_filename = '../data/KSTAR_detach_ctrl_Psi_data.h5'
efittree = 'EFIT01'

xpoint_margin = 0.2  # m
core_max_psin = 1
near_axis_psin = 0.2
show_figures = False
geometry = 'lower_null'

shots = [35853, 35854, 35857, 36161, 36162]

def get_psin_rz(shot):
    with h5py.File(efit_filename, 'r') as h5:
        psirz = h5[f'{shot}'][efittree]['PSIRZ']['data'][:]
        psia = h5[f'{shot}'][efittree]['SSIMAG']['data'][:]
        psib = h5[f'{shot}'][efittree]['SSIBRY']['data'][:]
        t = h5[f'{shot}'][efittree]['PSIRZ']['dim2'][:] * 1e-3
        r = h5[f'{shot}'][efittree]['PSIRZ']['dim0'][:]
        z = h5[f'{shot}'][efittree]['PSIRZ']['dim1'][:]
        zx1 = h5[f'{shot}'][efittree]['ZXPT1']['data'][:]
        rx2 = h5[f'{shot}'][efittree]['RXPT2']['data'][:]
        zx2 = h5[f'{shot}'][efittree]['ZXPT2']['data'][:]
        atime = h5[f'{shot}'][efittree]['RXPT1']['dim0'][:] * 1e-3
        rx1 = scipy.interpolate.interp1d(atime, h5[f'{shot}'][efittree]['RXPT1']['data'][:], bounds_error=False, fill_value=np.nan)(t)
        zx1 = scipy.interpolate.interp1d(atime, h5[f'{shot}'][efittree]['ZXPT1']['data'][:], bounds_error=False, fill_value=np.nan)(t)
        rx2 = scipy.interpolate.interp1d(atime, h5[f'{shot}'][efittree]['RXPT2']['data'][:], bounds_error=False, fill_value=np.nan)(t)
        zx2 = scipy.interpolate.interp1d(atime, h5[f'{shot}'][efittree]['ZXPT2']['data'][:], bounds_error=False, fill_value=np.nan)(t)

    psin = (psirz - psia[:, np.newaxis, np.newaxis])/((psib - psia)[:, np.newaxis, np.newaxis])

    return t, r, z, psin, rx1, zx1, rx2, zx2

def get_prad(shot):
    filename = f'../data/shot-{shot:06d}_data.mat'
    # Doesn't work because these files aren't v7.3
    # with h5py.File(filename, 'r') as mat:
    #     t = mat['time_Sec']
    #     power = mat['recon_MW']
    # Doesn't work because one array is 3d
    # data = mat4py.loadmat(filename)
    # t = data['time_Sec']
    # power = data['recon_MW']
    data = scipy.io.loadmat(filename)
    t = data['time_Sec'][0]
    power = data['recon_MW']
    prad_tot = data['Ptot_MW'][0]
    r = np.linspace(1.2555, 1.2555+0.011*99, 100)
    z = np.linspace(-1.4355, -1.4355+0.029*99, 100)
    return t, r, z, power, prad_tot

def analyze(shot, save=False):
    t, r, z, psin, rx1, zx1, rx2, zx2 = get_psin_rz(shot)
    tp, rp, zp, power_raw, prad_tot_official = get_prad(shot)
    t3d = t[:, np.newaxis, np.newaxis] + psin * 0
    z3d = z[np.newaxis, :, np.newaxis] + psin * 0
    r3d = r[np.newaxis, np.newaxis, :] + psin * 0
    power = scipy.interpolate.RegularGridInterpolator((tp, zp, rp), power_raw, bounds_error=False, fill_value=0.0)(np.array([t3d, z3d, r3d]).T).T
    #power[power < 0] = 0
    zxl = np.zeros(len(t))
    rxl = np.zeros(len(t))
    zxu = np.zeros(len(t))
    rxu = np.zeros(len(t))
    zxl[:] = np.nan
    zxu[:] = np.nan
    rxl[:] = np.nan
    rxu[:] = np.nan
    w1L = (-9 < zx1) & (zx1 < 0)
    w2L = (-9 < zx2) & (zx2 < 0)
    w1U = zx1 > 0
    w2U = zx2 > 0
    zxl[w1L] = zx1[w1L]
    rxl[w1L] = rx1[w1L]
    zxl[w2L] = zx2[w2L]
    rxl[w2L] = rx2[w2L]
    zxu[w1U] = zx1[w1U]
    rxu[w1U] = rx1[w1U]
    zxu[w2U] = zx2[w2U]
    rxu[w2U] = rx2[w2U]

    with h5py.File(efit_filename, 'r') as h5:
        lim = h5[f'{shot}'][efittree]['LIM']['data'][:]
        limx = lim[:, 0]
        limy = lim[:, 1]
    lim_polygon = shapely.geometry.Polygon(lim)
    within_lim2d = np.zeros((len(z), len(r)), bool)
    for ri, rr in enumerate(r):
        for zi, zz in enumerate(z):
            within_lim2d[zi, ri] = lim_polygon.contains(shapely.geometry.Point(rr, zz))

    within_lim = within_lim2d[np.newaxis, :, :] + (power * 0).astype(bool)
    power[~within_lim] = np.nan

    incore = (psin < core_max_psin) & within_lim
    indivl = (z3d < (zxl+xpoint_margin)[:, np.newaxis, np.newaxis]) & within_lim
    indivu = (z3d > (zxu-xpoint_margin)[:, np.newaxis, np.newaxis]) & within_lim
    nearaxis = (psin < near_axis_psin) & incore
    if geometry == 'lower_null':
        incore &= (z3d > (zxl+xpoint_margin)[:, np.newaxis, np.newaxis])
        indivu *= False
    elif geometry == 'upper_null':
        incore &= (z3d < (zxu-xpoint_margin)[:, np.newaxis, np.newaxis])
        indivu *= False
    else: # elif geometry == 'double_null':
        incore &= (z3d > (zxl+xpoint_margin)[:, np.newaxis, np.newaxis]) & (z3d < (zxu-xpoint_margin)[:, np.newaxis, np.newaxis])
    
    dr = r[1] - r[0]
    dz = z[1] - z[0]
    d = dr * dz * 2 * np.pi * r3d
    prad_core = np.nansum(d * power * incore, axis=(1,2))
    prad_axis = np.nansum(d * power * nearaxis, axis=(1,2))
    prad_divl = np.nansum(d * power * indivl, axis=(1,2))
    prad_divu = np.nansum(d * power * indivu, axis=(1,2))
    prad_tot = np.nansum(d * power, axis=(1,2))
    prad_sol = prad_tot - prad_core - prad_divl - prad_divu
    core_vol = np.sum(d * incore, axis=(1,2))
    divl_vol = np.sum(d * indivl, axis=(1,2))
    axis_vol = np.sum(d * nearaxis, axis=(1,2))

    if show_figures:
        tselect = 12
        i = np.argmin(abs(t-tselect))
        ip = np.argmin(abs(tp-tselect))
        power2d = power[i, :, :]
        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
        ax = axs[0, 0]
        ax.pcolormesh(r, z, power2d)
        ax.set_aspect('equal')
        ax = axs[1, 0]
        ax.pcolormesh(rp, zp, power_raw[ip])
        ax.set_aspect('equal')
        ax = axs[0, 1]
        power_div = copy.deepcopy(power)
        power_div[~indivl] = np.nan
        ax.pcolormesh(r, z, power_div[i, :, :])
        ax.set_aspect('equal')
        ax = axs[1, 1]
        power_core = copy.deepcopy(power)
        power_core[~incore] = np.nan
        ax.pcolormesh(r, z, power_core[i, :, :])
        core_power_slice = np.nansum(power_core[i, :, :])
        divl_power_slice = np.nansum(power_div[i, :, :])
        print(f'{core_power_slice=}, {divl_power_slice=}, {core_vol[i]=}, {divl_vol[i]=}, {axis_vol[i]=}')
        ax.set_aspect('equal')
        for ax in axs.flatten():
            ax.axhline(zxl[i]+xpoint_margin, color='white', ls='--')
            ax.axhline(zxu[i]-xpoint_margin, color='white', ls='--')
            ax.contour(r, z, psin[i, :, :], levels=[near_axis_psin, core_max_psin], colors=['white'])
            ax.plot(limx, limy, color='white')


    if show_figures:
        fig2, axs = plt.subplots(2)
        ax = axs[0]
        ax.plot(t, zx1, label='1')
        ax.plot(t, zx2, label='2')
        ax.plot(t, zxl, label='L', ls='--')
        ax.plot(t, zxu, label='U', ls='--')
        ax.legend()

        ax = axs[1]
        ax.plot(t, prad_core, label='core')
        ax.plot(t, prad_divl, label='divl')
        ax.plot(t, prad_divu, label='divu')
        ax.plot(t, prad_tot, label='total')
        ax.plot(tp, prad_tot_official, label='official total', color='k', ls='--')
        ax.plot(t, prad_sol, label='sol')
        ax.plot(t, prad_axis, label='near axis')
        ax.legend()

    if save:
        h5out_filename = f'../data/kstar_prad_analysis_{shot}.h5'
        with h5py.File(h5out_filename, 'w') as f:
            f.create_dataset('shot', data=shot)
            f.create_dataset('t', data=t)
            f.create_dataset('prad_tot', data=prad_tot)
            data_core = f.create_dataset('prad_core', data=prad_core)
            data_core.attrs['core_max_psin'] = core_max_psin
            data_core.attrs['xpoint_margin'] = xpoint_margin
            data_divl = f.create_dataset('prad_divl', data=prad_divl)
            data_divl.attrs['xpoint_margin'] = xpoint_margin
            data_divu = f.create_dataset('prad_divu', data=prad_divu)
            data_divu.attrs['xpoint_margin'] = xpoint_margin
            f.create_dataset('prad_sol', data=prad_sol)
            data_axis = f.create_dataset('prad_axis', data=prad_axis)
            data_axis.attrs['near_axis_psin'] = near_axis_psin
        print(f'Wrote {h5out_filename}')

    return t, prad_core, prad_divl, prad_divu



if show_figures:
    analyze(36161)
    plt.show()
else:
    for shot in shots:
        analyze(shot, save=True)