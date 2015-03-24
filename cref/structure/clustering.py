

def cluster_torsion_angles(torsion_angles):
    phi = torsion_angles['phi'][0]
    psi = torsion_angles['psi'][0]
    return {
        'H': (phi, psi),
        'E': (phi, psi),
        'C': (phi, psi),
    }
