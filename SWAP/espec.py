#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams
from scipy.io import readsav
from sys import argv
from numpy.ma import masked_array
from mom_trans_pluto import mom_trans

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

def makeespec(name):
    data = readsav(name+".sav")
    try:
        return _espec_v2(data, identify=_ident_any)
    except IndexError:
        print("Reading version 2 failed. Falling back on version 1.")
        return _espec_v1(data)

if __name__ == "__main__":
    fig, ax = makeespec(argv[1])

    mtrans = mom_trans(B0=0.08e-9)
    mtrans.plot_energies(ax,
            [mtrans.mp, mtrans.mpu, 2*mtrans.mp], 
            [{'color':'blue'},{'color':'red'},{'color':'green'}])

    ax.set_title("Synthetic SWAP Energy Spectrogram")

    fig.savefig(argv[1]+".png", format='png',bbox_inches='tight')
    fig.clear()
