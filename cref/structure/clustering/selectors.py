import logging

import numpy as np
from cref.structure.secondary import closest_ss

logger = logging.getLogger('CReF')


def top_scoring(model, scores):
    max_score = 0
    angles = []

    for k in range(model.n_clusters):
        my_members = model.labels_ == k
        my_scores = scores[my_members].tolist()
        score = sum(sorted(my_scores, reverse=True)[:5])
        if score > max_score:
            max_score = score
            angles = model.cluster_centers_[k]
    return angles


def secondary_structure(model, scores, structures, ss, ss_map):
    torsion_angles = secondary_structure_angles(
        model, structures, scores, ss_map)
    angles = select_cluster(torsion_angles, ss)
    return angles


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


def select_cluster(angles, ss):
    fallback_options = ('-', 'T', 'S', 'H', 'I', 'E', 'B')
    if ss in angles:
        return angles[ss]
    else:
        for fallback_ss in closest_ss(ss) + fallback_options:
            if fallback_ss in angles:
                logger.info(
                    'Could not find angles for {}, using angles for {}'.format(
                        ss, fallback_ss))
                return angles[fallback_ss]
