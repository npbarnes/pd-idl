#!/usr/bin/python
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams, gridspec
import colormaps as cmaps
plt.register_cmap(name='viridis', cmap=cmaps.viridis)
plt.register_cmap(name='plasma', cmap=cmaps.plasma)
from scipy.io import readsav
from sys import argv
from numpy.ma import masked_array
from numpy import logical_and
from mom_trans_pluto import mom_trans
import numpy as np
import spiceypy as sp
import NH_tools
from HybridReader2 import HybridReader2 as hr

def _espec_v1(data):
    x = data['xpos']
    ebins = data['ebins']
    cnt = data['cnt_arr']

    plt.pcolormesh(x,y,cnt,  norm=LogNorm())
    plt.yscale('log')
    plt.ylabel('Energy per charge ($\frac{eV}{q}$)')
    plt.xlabel('X position ($R_p$)')
    plt.title(title)
    plt.gca().invert_xaxis()
    cb = plt.colorbar()
    cb.set_label('Counts')

    return plt.gcf(), plt.gca()

def _ident_greater(cnt_arr, light, heavy):
    cnt_light = masked_array(cnt_arr,mask=(light <= heavy))
    cnt_heavy = masked_array(cnt_arr,mask=(light >  heavy))
    return cnt_light, cnt_heavy

def _ident_any(cnt_arr, light, heavy):
    cnt_light = masked_array(cnt_arr,mask=(heavy != 0))
    cnt_heavy = masked_array(cnt_arr,mask=(heavy == 0))
    return cnt_light, cnt_heavy


def _espec_v2(data, identify=None):
    if identify is None:
        identify = _ident_greater

    x = data['xpos']
    ebins = data['ebins']
    light = data['light_arr']
    heavy = data['heavy_arr']
    cnt_arr = data['cnt_arr']


    cnt_light, cnt_heavy = identify(cnt_arr, light, heavy)

    f = plt.gcf()
    ax = plt.gca()

    norm = LogNorm()

    lhist = ax.pcolormesh(x, ebins, cnt_light, norm=norm, cmap='Blues')
    lcb = plt.colorbar(lhist, fraction=0.075, pad=-0.02, shrink=0.75)
    lcb.ax.minorticks_on()

    hhist = ax.pcolormesh(x, ebins, cnt_heavy, norm=norm, cmap='Reds')
    hcb = plt.colorbar(hhist, fraction=0.075, pad=0.025, shrink=0.75, format="")

    lcb.ax.set_title('COIN (Hz)                     ', fontdict={'fontsize':'small'})
    ax.set_yscale('log')
    ax.set_ylabel('Energy per charge ($eV/q$)')
    ax.set_xlabel('X position ($R_p$)')
    ax.invert_xaxis()

    return f, ax

def _espec_v3(data, options=None):
    if options is None:
        options= {'colorbars':['Blues','Greens','Reds']}
    x = data['xpos']
    ebins = data['ebins']
    spectrograms_by_species = data['spectrograms_by_species']

    H = spectrograms_by_species[:,:,0]
    He = spectrograms_by_species[:,:,1]
    CH4 = spectrograms_by_species[:,:,2]

    cnt_arr = H+He+CH4

    mH = masked_array(cnt_arr, mask=(~logical_and(H>He,H>CH4)))
    mHe = masked_array(cnt_arr, mask=(~logical_and(He>H,He>CH4)))
    mCH4 = masked_array(cnt_arr, mask=~(CH4>0))

    f = plt.gcf()
    ax = plt.gca()

    norm = LogNorm()

    Hhist = ax.pcolormesh(x, ebins, mH, norm=norm, cmap=options['colorbars'][0])
    Hcb = plt.colorbar(Hhist, fraction=0.075, pad=-0.035, shrink=0.75)
    Hcb.ax.minorticks_on()

    Hehist = ax.pcolormesh(x, ebins, mHe, norm=norm, cmap=options['colorbars'][1])
    Hecb = plt.colorbar(Hehist, fraction=0.075, pad=-0.033, shrink=0.75, format="")

    CH4hist = ax.pcolormesh(x, ebins, mCH4, norm=norm, cmap=options['colorbars'][2])
    CH4cb = plt.colorbar(CH4hist, fraction=0.075, pad=0.025, shrink=0.75, format="")

    Hecb.ax.set_title('COIN (Hz)', fontdict={'fontsize':'small'})
    ax.set_yscale('log')
    ax.set_ylabel('Energy per charge ($eV/q$)')
    ax.set_xlabel('X position ($R_p$)')
    ax.invert_xaxis()

    return f, ax

def _espec_v4(data, options=None):
    if options is None:
        options= {'colorbars':['Blues','Greens','Reds']}

    # Load data from IDL
    spectrograms_by_species = data['spectrograms_by_species']
    t = data['times']
    x = data['positions'][0,:]
    y = data['positions'][1,:]
    z = data['positions'][2,:]
    o = data['orientations']
    ebins = data['ebins']

    # Separate parts
    H = spectrograms_by_species[:,:,0]
    He = spectrograms_by_species[:,:,1]
    CH4 = spectrograms_by_species[:,:,2]

    cnt_arr = H+He+CH4

    # Build masked arrays for plotting the dominant species
    mH = masked_array(cnt_arr, mask=(~logical_and(H>He,H>CH4)))
    mHe = masked_array(cnt_arr, mask=(~logical_and(He>H,He>CH4)))
    mCH4 = masked_array(cnt_arr, mask=~logical_and(CH4>H,CH4>He))
    
    # Setup subplot and colorbar axes
#    fig = plt.figure(figsize = (rcParams['figure.figsize'][0], 2*rcParams['figure.figsize'][1]))
#    gs = gridspec.GridSpec(8,1)
#    ax_spec  = plt.subplot(gs[0,0])
#    ax_theta = plt.subplot(gs[1,0], sharex=ax_spec)
#    ax_phi   = plt.subplot(gs[2,0], sharex=ax_spec)
#    ax_spin  = plt.subplot(gs[3,0], sharex=ax_spec)
#
#    ax_xy    = plt.subplot(gs[4:6,0])
#    ax_xz    = plt.subplot(gs[6:8,0], sharex=ax_xy)


    fig, (ax_spec, ax_theta, ax_phi, ax_spin, ax_xy, ax_xz) = plt.subplots(nrows=6, sharex=True, 
            gridspec_kw={'height_ratios':[1,1,1,1,2,3]}, 
            figsize=(rcParams['figure.figsize'][0], 2.25*rcParams['figure.figsize'][1]))

    fig.subplots_adjust(right=0.8, hspace=0.05)

    spec_pos = ax_spec.get_position()
    cbar_CH4 = fig.add_axes([spec_pos.x1+.01, spec_pos.y0+.01, 0.1/3, spec_pos.height-.02])
    cbar_CH4_pos = cbar_CH4.get_position()
    cbar_He = fig.add_axes([cbar_CH4_pos.x1, spec_pos.y0+.01, 0.1/3, spec_pos.height-.02])
    cbar_He_pos = cbar_He.get_position()
    cbar_H = fig.add_axes([cbar_He_pos.x1, spec_pos.y0+.01, 0.1/3,spec_pos.height-.02])

    # Plot spectrograms and make colorbars
    Hhist = ax_spec.pcolormesh(t, ebins, mH, norm=LogNorm(), cmap=options['colorbars'][0],
            vmin=1e-3, vmax=1e4)
    Hcb = fig.colorbar(Hhist, cax=cbar_H)
    Hcb.ax.minorticks_on()

    Hehist = ax_spec.pcolormesh(t, ebins, mHe, norm=LogNorm(), cmap=options['colorbars'][1],
            vmin=1e-3, vmax=1e4)
    Hecb = fig.colorbar(Hehist, cax=cbar_He, format="")

    CH4hist = ax_spec.pcolormesh(t, ebins, mCH4, norm=LogNorm(), cmap=options['colorbars'][2],
            vmin=1e-3, vmax=1e4)
    CH4cb = fig.colorbar(CH4hist, cax=cbar_CH4, format="")

    # Plot orientations
    ax_theta.plot(t, o[0], marker='o', markersize=1.3, linestyle='None', color='black')
    ax_phi.plot(t,   o[1], marker='o', markersize=1.3, linestyle='None', color='black')
    ax_spin.plot(t,  o[2], marker='o', markersize=1.3, linestyle='None', color='black')

    ax_theta.set_ylim(-180,180)
    ax_phi.set_ylim(-180,180)
    ax_spin.set_ylim(-180,180)

    # Plot positions
    h = hr('/home/nathan/data/2017-Wed-Aug-30/pluto-1/databig','np')
    n = h.get_last_timestep()[-1]

    # Get grid spacing
    qx = h.para['qx']
    qy = h.para['qy']
    qzrange = h.para['qzrange']

    # Find the center index of the grid
    cx = h.para['nx']/2
    cy = h.para['ny']/2
    cz = h.para['zrange']/2

    # the offset of pluto from the center isn't always availible
    try:
        po = h.para['pluto_offset']
    except KeyError:
        print("Couldn't get pluto_offset. It has been assumed to be 30, but it probably isn't.")
        po = 30

    # Set constatnt for Pluto radius 
    Rp = 1187. # km

    # Shift grid so that Pluto lies at (0,0,0) and convert from km to Rp
    qx = (qx - qx[len(qx)/2 + po])/Rp
    qy = (qy - qy[len(qy)/2])/Rp
    qzrange = (qzrange - qzrange[len(qzrange)/2])/Rp

    qt = [NH_tools.time_at_pos(x*Rp) for x in qx]

    T,Y = np.meshgrid(qt,qy)
    ax_xy.pcolormesh(T,Y,n[:,:,cz].transpose(), cmap='viridis', norm=LogNorm())

    T,Z = np.meshgrid(qt,qzrange)
    ax_xz.pcolormesh(T,Z,n[:,cy,:].transpose(), cmap='viridis', norm=LogNorm())

    ax_xy.plot(t, y/Rp-np.max(qy), color='black', linewidth=2)
    ax_xz.plot(t, np.zeros(z.shape), color='black', linewidth=2)

    # Set title, labels and adjust figure
    ax_spec.set_xlim(np.min(t), np.max(t))
    ax_spec.set_ylim(5e1, 1e4)
    ax_spec.set_title("Synthetic SWAP Energy Spectrogram", y=1.4)
    Hecb.ax.set_title('COIN (Hz)', fontdict={'fontsize':'small'})
    ax_spec.set_yscale('log')
    ax_spec.set_ylabel('Energy/Q ($eV/q$)')
    ax_theta.set_ylabel('Sun Theta\nAngle (deg)')
    ax_phi.set_ylabel('Sun Phi\nAngle (deg)')
    ax_spin.set_ylabel('Spin\nAngle (deg)')
    ax_xz.set_xlabel('Time (UTC)')

    ax_xy.set_ylabel('Y ($R_p$)')
    ax_xz.set_ylabel('Z ($R_p$)')
    ax_xy.set_ylim(np.min(qy), np.max(qy))
    ax_xz.set_ylim(1.5*np.min(qy), 1.5*np.max(qy))

    ax_xz.set_xticklabels([sp.timout(et, 'HR:MN:SC') for et in ax_spec.get_xticks()], rotation=-20, horizontalalignment='left')

    ax_twin = ax_spec.twiny()
    ax_twin.set_xlabel('X ($R_p$)')
    ax_twin.set_xlim(*ax_spec.get_xlim())
    ax_twin.set_xticks(ax_spec.get_xticks())
    ax_twin.set_xticklabels(['{0:.2f}'.format(NH_tools.pos_at_time(et)/1187) for et in ax_spec.get_xticks()], y=1.05)

    # This breaks the convention of previous versions
    return fig, None


def makeespec(name, options=None):
    data = readsav(name+".sav")
    try:
        return _espec_v4(data, options)
    except KeyError:
        try:
            print("Reading version 4 failed. Falling back on version 3.")
            return _espec_v3(data, options)
        except KeyError:
            try:
                print("Reading version 3 failed. Falling back on version 2.")
                return _espec_v2(data, identify=_ident_any)
            except IndexError:
                print("Reading version 2 failed. Falling back on version 1.")
                return _espec_v1(data)

if __name__ == "__main__":
    fig, ax = makeespec(argv[1], {'colorbars':['Blues','Greens','Reds']})

    #mtrans = mom_trans(B0=0.3e-9)
    #mtrans.plot_energies(ax,
    #        [mtrans.mp, mtrans.mpu, 2*mtrans.mp], 
    #        [{'color':'blue'},{'color':'red'},{'color':'green'}])


    fig.savefig(argv[1]+".png", format='png', bbox_inches='tight')
    fig.clear()
