#!/usr/bin/python
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams
from scipy.io import readsav
from sys import argv
from numpy.ma import masked_array
from numpy import logical_and
from mom_trans_pluto import mom_trans
import numpy as np
import spiceypy as sp
import NH_tools

def _espec_v1(data):
    x = data['xpos']
    ebins = data['ebins']
    cnt = data['cnt_arr']

    plt.pcolormesh(x,y,cnt, norm=LogNorm())
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

    spectrograms_by_species = data['spectrograms_by_species']
    t = data['times']
    x = data['positions'][0,:]/1187
    o = data['orientations']
    ebins = data['ebins']

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

    Hhist = ax.pcolormesh(t, ebins, mH, norm=norm, cmap=options['colorbars'][0])
    Hcb = plt.colorbar(Hhist, fraction=0.075, pad=-0.035, shrink=0.75)
    Hcb.ax.minorticks_on()

    Hehist = ax.pcolormesh(t, ebins, mHe, norm=norm, cmap=options['colorbars'][1])
    Hecb = plt.colorbar(Hehist, fraction=0.075, pad=-0.033, shrink=0.75, format="")

    CH4hist = ax.pcolormesh(t, ebins, mCH4, norm=norm, cmap=options['colorbars'][2])
    CH4cb = plt.colorbar(CH4hist, fraction=0.075, pad=0.025, shrink=0.75, format="")


    ax.set_title("Synthetic SWAP Energy Spectrogram", y=1.11)
    Hecb.ax.set_title('COIN (Hz)', fontdict={'fontsize':'small'})
    ax.set_yscale('log')
    ax.set_ylabel('Energy per charge ($eV/q$)')
    ax.set_xlabel('Time (UTC)')

    ax.set_xticklabels([sp.timout(et, 'HR:MN:SC') for et in ax.get_xticks()], rotation=30)

    ax2 = ax.twiny()
    ax2.set_xlabel('X ($R_p$)')
    ax2.set_xlim(*ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(['{0:.2f}'.format(NH_tools.pos_at_time(et)/1187) for et in ax.get_xticks()], y=1.01)

    return f, ax


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
