import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing
import matplotlib.pyplot as plt
from scipy import sparse

from cref.structure.plot import ramachandran_surface
from cref.structure.secondary import ss_eight_to_three, closest_ss


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


def secondary_structure_angles(model, structures, all_scores, ss_eight=True):
    angles = {}
    scores = {}
    for k in range(model.n_clusters):
        my_members = model.labels_ == k
        my_structures = structures[my_members].tolist()
        my_scores = all_scores[my_members].tolist()
        ss = max(set(my_structures), key=my_structures.count)
        score = sum(sorted(my_scores, reverse=True)[:5])
        if not ss_eight:
            ss = ss_eight_to_three(ss)
        if (ss not in angles) or (scores[ss] < score):
            angles[ss] = model.cluster_centers_[k]
            scores[ss] = score
    return angles


def select_cluster(angles, ss, select):
    if ss in angles:
        return angles[ss]
    else:
        print('Could not find angles for ss ' + ss)
        for fallback_ss in closest_ss(ss) + ('C', 'H', 'E'):
            if fallback_ss in angles:
                print('Using angles for {} instead'.format(fallback_ss))
                return angles[fallback_ss]


def cluster_torsion_angles(blast_structures, ss, n_clusters=6,
                           select="ss", ss_eight=True):
    phi = blast_structures['phi']
    psi = blast_structures['psi']
    structures = blast_structures['central_ss']

    # Encode each secondary structure as the ascii value for the character
    enc = OneHotEncoder()
    ord_structures = [[ord(c)] for c in structures]
    enc.fit(ord_structures)
    encoded_structures = enc.transform(ord_structures)
    scores = blast_structures['score']

    phi_scaler = preprocessing.StandardScaler().fit(phi)
    psi_scaler = preprocessing.StandardScaler().fit(psi)

    X = np.vstack((zip(phi_scaler.transform(phi), psi_scaler.transform(psi))))
    X = sparse.hstack((X, encoded_structures))

    model = KMeans(init='k-means++', n_clusters=n_clusters)
    model.fit(X)
    torsion_angles = secondary_structure_angles(
        model, structures, scores, ss_eight)

    # plot_clusters(model, X, blast_structures['fragment'][0])
    angles = select_cluster(torsion_angles, ss, select)
    angles = (
        phi_scaler.inverse_transform(angles[0]),
        psi_scaler.inverse_transform(angles[1])
    )
    return angles
