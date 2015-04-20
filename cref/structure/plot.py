from collections import OrderedDict

import matplotlib.pyplot as plt
import pandas


_ramachandran_densities = pandas.read_csv(
    'data/rama500-general.data',
    skiprows=6,
    delimiter=' ',
    names=['phi', 'psi', 'value']
)

"""
DSSP output:
    H = α-helix
    B = residue in isolated β-bridge
    E = extended strand, participates in β ladder
    G = 3-helix (310 helix)
    I = 5 helix (π-helix)
    T = hydrogen bonded turn
    S = bend

Colors extracted from rcsb.org.
"""

DSSP_to_color = {
    'H': '#ED6161',
    'B': '#CCA200',
    'E': '#FFFB00',
    'G': '#FFC2C2',
    'I': '#900000',
    'T': '#990099',
    'S': '#0000FF',
    '-': 'black',
}


def ramachandran_surface():
    """
    Plot density surface for generic ramachandran
    """
    fontsize = 18
    ticks = [-180, -90, 0, 90, 180]
    plt.contourf(
        list(OrderedDict.fromkeys(_ramachandran_densities['phi'])),
        list(OrderedDict.fromkeys(_ramachandran_densities['psi'])),
        _ramachandran_densities['value'].reshape(180, 180).T,
        levels=[0, 0.0005, 0.02, 1],
        colors=['#FFFFFF', '#B3E8FF', '#7FD9FF']
    )
    plt.xlabel('$\phi$', fontsize=fontsize)
    plt.ylabel('$\psi$', fontsize=fontsize)
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.tick_params(direction="out")
    plt.margins(0.05)
    ax = plt.axes()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def ramachandran(torsion_angles, fragment, central):
    """
    Plot ramachandran of a set of torsion angles for a given fragment

    :param torsion_angles: Dictionary with torsion angles phi and psi
    :param fragment: Fragment identifier, used for displaying purposes
    """
    ramachandran_surface()
    plt.title("Ramachandran plot for " + fragment)
    plt.scatter(
        x=torsion_angles['phi'],
        y=torsion_angles['psi'],
        s=[1.05 ** x for x in torsion_angles['identity']],
        c=[DSSP_to_color[ss] for ss in torsion_angles['central_ss']],
        cmap="seismic",
        marker='o',
        alpha=0.5,
    )
    plt.show()
