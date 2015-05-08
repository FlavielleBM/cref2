import logging

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing
import matplotlib.pyplot as plt
from scipy import sparse

from cref.structure.plot import ramachandran_surface
from cref.structure.secondary import closest_ss

logger = logging.getLogger('CReF')


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


def secondary_structure_angles(model, structures, all_scores, ss_map):
    angles = {}
    scores = {}
    for k in range(model.n_clusters):
        # tie breaker on the top 5 score
        my_members = model.labels_ == k
        my_scores = all_scores[my_members].tolist()
        score = sum(sorted(my_scores, reverse=True)[:5])

        ss = chr(ss_map[np.argmax(model.cluster_centers_[k][2:])])

        if (ss not in angles) or (scores[ss] < score):
            angles[ss] = model.cluster_centers_[k]
            scores[ss] = score
    return angles


def select_cluster(angles, ss, select):
    if ss in angles:
        return angles[ss]
    else:
        for fallback_ss in closest_ss(ss) + ('C', 'H', 'E'):
            if fallback_ss in angles:
                logger.info(
                    'Could not find angles for {}, using angles for {}'.format(
                        ss, fallback_ss))
                return angles[fallback_ss]


def cluster_torsion_angles(blast_structures, ss, n_clusters=8, select="ss"):
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
    torsion_angles = secondary_structure_angles(
        model, structures, scores, ss_map)

    # plot_clusters(model, X, blast_structures['fragment'][0])
    angles = select_cluster(torsion_angles, ss, select)
    logger.info('Selected cluster: {} {}'.format(
        ss, angles[2:], chr(enc.active_features_[np.argmax(angles[2:])])))
    angles = (
        phi_scaler.inverse_transform(angles[0]),
        psi_scaler.inverse_transform(angles[1])
    )
    return angles
