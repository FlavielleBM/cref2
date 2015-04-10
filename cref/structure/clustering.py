

def cluster_torsion_angles(blast_structures):
    phi = blast_structures['phi'][0]
    psi = blast_structures['psi'][0]
    return {
        'H': (phi, psi),
        'E': (phi, psi),
        'C': (phi, psi),
    }
