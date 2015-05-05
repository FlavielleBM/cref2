import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

from cref.structure.plot import ramachandran_surface
from cref.structure.secondary import dssp_to_porter


def plot_clusters(model, X, fragment):
    ramachandran_surface()
    colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan']
    for k, col in zip(range(model.n_clusters), colors):
        my_members = model.labels_ == k
        cluster_center = model.cluster_centers_[k]
        plt.plot(X[my_members, 0], X[my_members, 1], 'w',
                 markerfacecolor=col, marker='*', markersize=5)
        plt.plot(cluster_center[0], cluster_center[1], 'o',
                 markerfacecolor=col, markeredgecolor='k', markersize=8)
    plt.title('KMeans for fragment ' + fragment)
    plt.show()


def secondary_structure_angles(model, structures):
    angles = {}
    for k in range(model.n_clusters):
        my_members = model.labels_ == k
        my_structures = structures[my_members].tolist()
        ss = max(set(my_structures), key=my_structures.count)
        porter_ss = dssp_to_porter(ss)
        if porter_ss not in angles:
            angles[porter_ss] = model.cluster_centers_[k]
    return angles


def cluster_torsion_angles(blast_structures, ss, n_clusters=6):
    phi = blast_structures['phi']
    psi = blast_structures['psi']
    structures = blast_structures['central_ss']
    # score = blast_structures['score']

    X = np.vstack((zip(phi, psi)))

    model = KMeans(init='k-means++', n_clusters=n_clusters, n_init=10)
    model.fit(X)

    # plot_clusters(model, X, blast_structures['fragment'][0])
    angles = secondary_structure_angles(model, structures)
    if ss in angles:
        return angles[ss]
    else:
        print('Could not find angles for ss ' + ss)
        return angles.values()[0]
