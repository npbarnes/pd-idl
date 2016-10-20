#!/usr/bin/python3
from espec import makeespec

for name in ["espec-"+str(i+1) for i in range(8)]:
    fig = makeespec(name, "All particles")

    fig.savefig(name+".png", format='png',bbox_inches='tight')
    fig.clear()

