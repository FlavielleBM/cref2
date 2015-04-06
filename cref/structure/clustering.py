

def cluster_torsion_angles(torsion_angles):
    phi = torsion_angles[0][1]
    psi = torsion_angles[0][2]
    return {
        'H': (phi, psi),
        'E': (phi, psi),
        'C': (phi, psi),
    }
