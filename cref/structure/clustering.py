import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

from cref.structure.plot import ramachandran_surface


def plot_clusters(model, X, fragment):
    ramachandran_surface()
    colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan']
    for k, col in zip(range(model.n_clusters), colors):
        my_members = model.labels_ == k
        cluster_center = model.cluster_centers_[k]
        plt.plot(X[my_members, 0], X[my_members, 1], 'w',
                 markerfacecolor=col, marker='*', markersize=5)
        plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=8)
    plt.title('KMeans for fragment ' + fragment)
    plt.show()


def cluster_torsion_angles(blast_structures, n_clusters=6):
    phi = blast_structures['phi']
    psi = blast_structures['psi']
    X = np.vstack((zip(phi, psi)))

    model = KMeans(init='k-means++', n_clusters=6, n_init=10)
    model.fit(X)
    plot_clusters(model, X, blast_structures['fragment'][0])

    angles = model.cluster_centers_[0]
    return {
        'H': angles,
        'E': angles,
        'C': angles,
    }
