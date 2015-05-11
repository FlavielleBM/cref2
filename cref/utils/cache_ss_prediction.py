from cref.structure.secondary import SecondaryStructurePredictor


def cache_ss_prediction(sspred, order):
    curs = sspred.pdb_db.conn.cursor().execute(
        'select distinct sequence from pdb_ss order by pdb ' + order)
    seqs = curs.fetchall()
    for seq in seqs:
        sspred.sspro(seq[0])

if __name__ == '__main__':
    import sys
    sspred = SecondaryStructurePredictor(sys.argv[2])
    cache_ss_prediction(sspred, sys.argv[1])
