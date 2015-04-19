blastdb:
	wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz -O tests/blastdb/pdb_seqres.txt.gz
	cd tests/blastdb && gunzip -f pdb_seqres.txt.gz && \
	makeblastdb -in pdb_seqres.txt -out pdbseqres -dbtype prot	

torsions:
	cc -o cref/structure/torsions cref/structure/torsions.c -lm

peptide:
	git clone https://github.com/mchelem/peptide

python_packages:
	virtualenv --python=/usr/bin/python3 venv
	bash -i -c "source venv/bin/activate && pip install -r requirements.txt"

packages:
	sudo apt-get install -y tcl-dev tk-dev
	sudo apt-get install -y liblapack-dev gcc gfortran
	sudo apt-get install -y ncbi-blast+
	sudo apt-get install -y python3 python3-dev python-virtualenv

ss:
	wget http://www.rcsb.org/pdb/files/ss.txt.gz -O data/ss.txt.gz
	cd data && gunzip -f ss.txt.gz
	bash -i -c "source venv/bin/activate && python -m cref.utils.import_pdb_ss data/ss.txt data/ss.db"

install: packages python_packages peptide torsions blastdb ss
