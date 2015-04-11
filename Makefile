update_blastdb:
	wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz -O tests/blastdb/pdb_seqres.txt.gz
	cd tests/blastdb && gunzip pdb_seqres.txt.gz && \
	makeblastdb -in pdb_seqres.txt -out pdbseqres -dbtype prot	

create_blastdb:
	cd tests/blastdb && gunzip pdb_seqres.txt.gz && \
	makeblastdb -in pdb_seqres.txt -out pdbseqres -dbtype prot	

torsions:
	cc -o cref/structure/torsions cref/structure/torsions.c -lm

peptide:
	git clone git@github.com:mchelem/peptide.git

