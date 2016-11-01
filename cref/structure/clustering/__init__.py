import logging
import os

from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing

import matplotlib.pyplot as plt
from pylab import cm

from scipy import sparse
import numpy as np


from cref.structure.clustering import selectors
from cref.structure.plot import ramachandran_surface

logger = logging.getLogger('CReF')


def plot_clusters(model, SX, angles, phi_scaler, psi_scaler, fragment,
                  output_writer=None, output_dir=None):
    ramachandran_surface()
    X = SX.todense()
    # Start from 10 to avoid blacks
    colors = np.linspace(10, 255, model.n_clusters)
    for k, color in zip(range(model.n_clusters), colors):
        col = cm.spectral(int(color))
        my_members = model.labels_ == k
        # cluster_center = model.cluster_centers_[k]
        plt.scatter(
            x=phi_scaler.inverse_transform(X[my_members, 0]),
            y=psi_scaler.inverse_transform(X[my_members, 1]),
            c=col,
            s=30,
            marker='o',
            alpha=0.5
        )
        plt.title('Clusters for {} (inertia: {:.2})'.format(
            fragment, model.inertia_))
    if output_writer:
        output_writer.savefig(dpi=150)
    if output_dir:
        plt.savefig(
            os.path.join(output_dir, 'clustering', fragment + '.cluster.svg'),
            format='svg', dpi=300,
        )
    plt.close()


def cluster_torsion_angles(blast_structures, ss, n_clusters=8,
                           selector="ss", name="cluster_plot",
                           output_writer=None, output_dir=None):
    phi = blast_structures['phi'].values.reshape(-1, 1)
    psi = blast_structures['psi'].values.reshape(-1, 1)
    structures = blast_structures['central_ss']

    # Encode each secondary structure as the ascii value for the character
    enc = OneHotEncoder()
    ord_structures = [[ord(c)] for c in structures]
    enc.fit(ord_structures)
    encoded_structures = enc.transform(ord_structures)
    ss_map = enc.active_features_

    scores = blast_structures['score']

    # Scale phi and psi to a standard normal
    phi_scaler = preprocessing.StandardScaler().fit(phi)
    psi_scaler = preprocessing.StandardScaler().fit(psi)

    X = np.hstack((phi_scaler.transform(phi), psi_scaler.transform(psi)))
    X = sparse.hstack((X, encoded_structures))
    print("SHAPEFIM", X.shape)

    model = KMeans(init='k-means++', n_clusters=n_clusters, random_state=1)
    model.fit(X)
    if selector == "score":
        angles = selectors.top_scoring(model, scores)
    else:
        angles = selectors.secondary_structure(
            model, scores, structures, ss, ss_map)
    inertia = model.inertia_
    logger.info('Selected cluster: {} {}'.format(
        ss, angles[2:], chr(enc.active_features_[np.argmax(angles[2:])])))
    plot_clusters(
        model, X, angles,
        phi_scaler, psi_scaler,
        name, output_writer, output_dir,
    )
    angles = (
        phi_scaler.inverse_transform(angles[0].reshape(-1, 1)),
        psi_scaler.inverse_transform(angles[1].reshape(-1, 1))
    )
    return angles, inertia
