import logging

from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing
import matplotlib.pyplot as plt
from scipy import sparse
import numpy as np

from cref.structure.clustering import selectors
from cref.structure.plot import ramachandran_surface

logger = logging.getLogger('CReF')


def plot_clusters(model, SX, angles, phi_scaler, psi_scaler, fragment):
    plt.figure()
    ramachandran_surface()
    colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan']
    X = SX.todense()
    for k, col in zip(range(model.n_clusters), colors):
        my_members = model.labels_ == k
        cluster_center = model.cluster_centers_[k]
        plt.plot(
            phi_scaler.inverse_transform(X[my_members, 0]),
            psi_scaler.inverse_transform(X[my_members, 1]),
            'w', markerfacecolor=col, marker='*', markersize=5
        )
        if angles[0] == cluster_center[0]:
            plt.plot(
                phi_scaler.inverse_transform(cluster_center[0]),
                psi_scaler.inverse_transform(cluster_center[1]),
                'D', markerfacecolor=col, markeredgecolor='k', markersize=8
            )
        else:
            plt.plot(
                phi_scaler.inverse_transform(cluster_center[0]),
                psi_scaler.inverse_transform(cluster_center[1]),
                'o', markerfacecolor=col, markeredgecolor='k', markersize=8
            )
    plt.title('KMeans for fragment ' + fragment)
    plt.savefig('predictions/tmp/{}_clustering.png'.format(fragment))
    plt.close()


def cluster_torsion_angles(blast_structures, ss, n_clusters=8, selector="ss"):
    phi = blast_structures['phi']
    psi = blast_structures['psi']
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

    X = np.vstack((zip(phi_scaler.transform(phi), psi_scaler.transform(psi))))
    X = sparse.hstack((X, encoded_structures))

    model = KMeans(init='k-means++', n_clusters=n_clusters)
    model.fit(X)
    if selector == "score":
        angles = selectors.top_scoring(model, scores)
    else:
        angles = selectors.secondary_structure(
            model, scores, structures, ss, ss_map)
    inertia = model.inertia_
    logger.info('Selected cluster: {} {}'.format(
        ss, angles[2:], chr(enc.active_features_[np.argmax(angles[2:])])))
    plot_clusters(model, X, angles, phi_scaler, psi_scaler,
                  blast_structures['fragment'][0])
    angles = (
        phi_scaler.inverse_transform(angles[0]),
        psi_scaler.inverse_transform(angles[1])
    )
    return angles, inertia
