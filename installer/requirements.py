 #!/usr/bin/python3
######################################################
#####################REQUIREMENTS#####################
######################################################
pythonVer  = 'python3.6'

systemPkgs = [pythonVer, '%s-dev' % (pythonVer), 'python3-tk', 'virtualenv', 'git', 
            'tcl8.5-dev', 'tk8.5-dev', 'liblapack-dev', 'gcc', 'gfortran', 'ncbi-blast+', 'pymol']

pipPkgs    = ['pandas', 'requests', 'biopython==1.70', 'scikit-learn==0.19', 'scipy',
            'matplotlib==1.5.3', 'openpyxl', 'git+https://github.com/mchelem/porter_paleale']
######################################################
######################################################
######################################################