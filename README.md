CReF
=====

Implementation of CReF: a Central Residue Fragment-based method to predict approximate 3-D structures of polypeptides by mining the Protein Data Bank (PDB).

## Installation

```make install```

The command downloads and builds dependencies. For details, read the Makefile.

## Usage

### List available arguments:

```./cref.sh --help```


### Predict a protein from the pdb and display it using pymol:

```./cref.sh --pdb 1zdd --pymol```

By default the resulting PDB is saved to the predictions folder and named predicted_structure.pdb

### Predict a protein from an input sequence and output to a given directory:

```./cref.sh --sequence DFYFNAI --output predictions/myprotein/```

### Predict a protein from a fasta file:


```/cref.sh --fasta 1zdd.fasta```

## Configuration

The [config.py](https://github.com/mchelem/cref2/blob/master/config.py) file contains the default configuration. You can change parameters and specify it on the command line:

```./cref.sh --config config --pdb 1zdd --pymol```

Note you MUST NOT include the .py in the command line.

The BLAST parameters are used with the ncbi-blast tool. See https://www.ncbi.nlm.nih.gov/books/NBK279675/ for more information.




