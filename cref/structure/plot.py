from collections import OrderedDict

import matplotlib.pyplot as plt
import pandas


_ramachandran_densities = pandas.read_csv(
    'data/rama500-general.data',
    skiprows=6,
    delimiter=' ',
    names=['phi', 'psi', 'value']
)

def ramachandran_surface():
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


def ramachandran(torsion_angles, fragment):
    """
    Plot ramachandran of a set of torsion angles for a given fragment

    :param torsion_angles: Dictionary with torsion angles phi and psi
    :param fragment: Fragment identifier, used for displaying purposes
    """
    ramachandran_surface()
    plt.title("Ramachandran plot for " + fragment)
    plt.scatter(
        torsion_angles['phi'], torsion_angles['psi'],  cmap="b", marker='.')
    plt.show()
