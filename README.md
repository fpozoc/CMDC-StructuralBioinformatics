# CMDC.py
## Developed by Fernando Pozo Ocampo & Marcos Camara Donoso
## PROJECT nº3: "CORRELATED MUTATIONS AND DISTANCE CORRELATIONS TO PREDICT AMINOACID INTERACTIONS"
## SUBJECTS: Structural Bioinformatics and Introduction to Python (Baldo Oliva, Javier García)
## March 12, 2017
======
report.pdf contains all of execution and software details. Check section 2. Methods: CMDC.py Step by Step in order to know more about the usage of CMDC.py.
======
For a fast execution, type in the command line on this directory:
$ chmod u+x run.sh
$ ./run.sh
After approximately 10 minutes, you will find in your directory:
$ ls
CMDCresults run.sh README.md CMDC.py distances.py mi.py extract_sequences.py 
$ cd CMDCresults
$ ls
CMDC_01Example_results CMDC_02Example_results CMDC_03Example_results 
$ cd CMDC_01Example_results
$ ls
my_scripts outputs pdb5cyt.ent plots std.sys

NOTE: run.sh have already executed in this directory and you can find the results in directory CMDCresults. If you want to execute it again, please, remove CMDCresults from this directory in order to get new results.
======
If you want to run this 3 examples independiently, type in the command line:
$ python3 CMDC.py -atom CB 5cyt -seqs 600 -gaps -msa clustalw
$ python3 CMDC.py 1rbb -gaps -msa muscle
$ python3 CMDC.py 3bp2 -a CA -seqs 50 -msa clustalo 

Type:
$ python CMDC.py -h 
If you want to know more about the optional arguments of the program.
