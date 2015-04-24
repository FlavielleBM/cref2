blastdb:
	mkdir -p data/blastdb
	wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz -O data/blastdb/pdb_seqres.txt.gz
	cd data/blastdb && gunzip -f pdb_seqres.txt.gz && \
	makeblastdb -in pdb_seqres.txt -out pdbseqres -dbtype prot	
	rm -rf data/pdb_seqres.txt

pdb:
	mkdir -p data/pdb_rsync
	mkdir -p data/pdb
	rsync -rlpt -v -z --delete --port=33444 \
		rsync.wwpdb.org::ftp_data/structures/divided/pdb/ data/pdb_rsync/
	cp -r data/pdb_rsync/* data/pdb/
	gunzip -r data/pdb/

torsions:
	cc -o cref/structure/torsions cref/structure/torsions.c -lm

peptide:
	git clone https://github.com/mchelem/peptide

python_packages:
	virtualenv --python=/usr/bin/python3 env
	bash -i -c "source env/bin/activate && pip install -r requirements.txt"

packages:
	sudo apt-get install -y tcl8.4-dev tk8.4-dev libpng3
	sudo apt-get install -y liblapack-dev gcc gfortran
	sudo apt-get install -y ncbi-blast+
	sudo apt-get install -y rabbitmq-server
	sudo apt-get install -y python3 python3-dev python3-tk python-virtualenv

ss:
	wget http://www.rcsb.org/pdb/files/ss.txt.gz -O data/ss.txt.gz
	cd data && gunzip -f ss.txt.gz
	bash -i -c "source env/bin/activate && python -m cref.utils.import_pdb_ss data/ss.txt data/ss.db"
	rm -rf data/ss.txt

install: packages python_packages peptide torsions blastdb ss
	mkdir predictions
