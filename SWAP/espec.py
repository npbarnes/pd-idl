#!/usr/bin/python3
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.io import readsav
from sys import argv

def makeespec(name, title):
    dat = readsav(name+".sav")
    x = dat['xpos']
    y = dat['ebins']
    cnt = dat['cnt_arr']

    plt.pcolormesh(x,y,cnt, norm=LogNorm())
    plt.yscale('log')
    plt.ylabel('Energy per charge (eV/q)')
    plt.xlabel('X position (R_p)')
    plt.title(title)
    plt.gca().invert_xaxis()
    cb = plt.colorbar()
    cb.set_label('Counts')

    return plt.gcf()

    plt.savefig(name+".png", format='png',bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    fig = makeespec(argv[1],argv[2])
    fig.savefig(argv[1]+".png", format='png',bbox_inches='tight')
    fig.clear()
