"""DOCUMENTATION FOR CMDC.py: "CORRELATION MUTATIONS AND DISTANCE CORRELATIONS TO PREDICT AMINOACIDS INTERACTION". \nWritten in Python3 by Marcos Camara Donoso & Fernando Pozo Ocampo. Project for the subjects Introduction to Python and Structural Bioinformatics. M.Sc. in Bionformatics for Health Sciences. Pompeu Fabra University.\n\n
extract_sequences.py standalone module.\n
This module runs a BLAST alignment provided by the NCBI against our protein FASTA sequence. The sequences are retrieved and a multiple sequence alignment (MSA) is performed. User can select into 4 different types of MSA (ClustalOmega, ClustalW, T-Coffee and Muscle).\n\n module sys, argparse, os, urllib and BioPython are required by this part of the program."""

import sys, argparse, os, urllib, ftplib
from distutils.spawn import find_executable
from Bio import Entrez, SeqIO # Return the data as a handle object.
                              # Standard input - output interface in Biopython
from Bio.Blast import NCBIWWW, NCBIXML # NCBIWWW for running BLAST online # NCBIXML es el modulo más sencillo y usado para parsear el output
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1  # Turns a three-letter code protein sequence into one with single letter codes
from Bio.PDB import PDBList
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline, TCoffeeCommandline, ClustalOmegaCommandline
from Bio.Application import ApplicationError


# CHECK THIS DOCUMENTATION
# http://biopython.org/DIST/docs/api/Bio.Blast.NCBIWWW-pysrc.html#
# http://biopython.org/DIST/docs/api/Bio.Blast.NCBIXML-pysrc.html#
def BLAST_online(query, blast, database, hits):
    sys.stdout.write("BLAST ALIGNMENT: query:{}, database:{}, type:{}.\n".format(query,database,blast))

    if os.path.isfile(query): # This 'if statement' checks if your fasta file is in the same directory of the program your are running
        record = SeqIO.read(query, format="fasta") # Reading input in fasta format
        result_handle = NCBIWWW.qblast(blast, database, record.format("fasta"), hitlist_size=hits) # Hitlist_size is number of hits to return (in blast module default is 50)
    else: # If you don't have a fasta file locally, check the PDB id you are introducing
        result_handle = NCBIWWW.qblast(blast, database, query, hitlist_size=hits)
    blast_record = NCBIXML.read(result_handle) # Save the output in a variable called blast_record. Format XML is needed in order to be parsed

    return blast_record

# Y para imprimir imprimir el xml...
#my_query = SeqIO.read("pdbidentifier.fasta", format="fasta")
#result_handle = NCBIWWW.qblast("blastp", "nr", my_query.seq)
#blast_printing_result = open("my_blast.xml", "w")
#blast_printing_result.write(result_handle.read())
#blast_result.close()
#result_handle.close()

# Forma clásica que añadimos aquí para el parseo de un blast xml (http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc93)
# Asumimos que en casos de fuertes correlación a nivel de aminoácidos es menos posible para reflejar para conser
def BLAST_filter(blast_handle, genus=True):
    """It analyse the result of BLAST alignment. You can filter your result by gender. By default, our option is not to filt."""
    sys.stdout.write("The ID you are try to align is {}.\n".format(blast_handle.query_id))

    length_query = blast_handle.query_length  # Checking lenght of the alignment in order to select the create a percentage_identity
    identifier_list = [] # list of identifiers
    gen_set = set()

    for alignment in blast_handle.alignments:  # for each BLAST alignment
        sys.stdout.write("ALIGNMENT:{}\n".format(alignment.title))

# You can get better results if filter by gender in BLAST results parse.
# https://www.ncbi.nlm.nih.gov/pubmed/11707606 Similarity of phylogenetic trees as indicator of protein-protein interaction. Valencia y Pazos lo describen en el 2001
        if genus: # Selecting filter by gender
            if "[" in alignment.title:# title attribute is tame of alignment, the genere is into always appears into [] after name of the protein and its identifiers
                specie = alignment.title.split('[')[-1].rstrip("]") # Splitting the name by genere and specie. Takes the second (-1=last) name of the list
                gen = specie.split()[0] # Saving the genere in gen. It is the first element of the list
# CUIDADO CON EL PARSEO DEL XML PARA SECUENCIAS EXTRAÍDAS DE BASES DE DATOS DIFERENTES A NR
                for hsp in alignment.hsps: # Iterating each protein hit
                    perc_ident = 100 * hsp.identities / length_query
# LEER AQUÍ ACLARACIÓN DE POR QUÉ HACEMOS ESTO ---> Quantification of the variation in percentage identity for protein sequence alignments. Raghava and Barton, 2006
                    if perc_ident > 30:
                        if gen not in gen_set:
                            seps = alignment.hit_id.split("|") # Splitting the 4 values of hit identifier
                            id_add = {seps[0]: seps[1],
                                      seps[2]: seps[3]}
                            identifier_list.append(id_add) # Creating the list of dictionaries
                            gen_set.add(gen)# Adding to a set diferent generes each time
                    else:
                        sys.stdout.write("FINISHING ALIGNMENT WITH HIGH SCORE...\n")
        else: # If you don't want to filter by genere
            seps = alignment.hit_id.split("|")
            identifier_list.append({seps[0]: seps[1], seps[2]: seps[3]}) # final list of identifiers
    else:
        sys.stdout.write("MY ALIGNMENT HAS ALREADY FINISHED.\n")

    return identifier_list

# Entrez te pide que un email y el nombre de la herramienta que estás usando para coger tu PDB
Entrez.tool = "CMDC_project.py"
Entrez.email = "fernando.pozo01@estudiant.upf.edu"

# http://biopython.org/DIST/docs/api/Bio.Entrez-module.html#efetch
def extract_sequence_from_BLASTids(identifier_seqs):
    sys.stdout.write("DOWNLOADING SEQUENCES\n")

    for identifier in identifier_seqs:
        sys.stdout.write("RETRIEVE {}.\n".format(identifier))

        efetch_filehandle = Entrez.efetch(db="protein", id=identifier, rettype="fasta", retmode="text") # Retrieves sentences of proteins in fasta format from a list of one or more primary IDs or from the user's environment
        yield SeqIO.read(efetch_filehandle, "fasta") # Generating all of fasta sentences I need


def filt_ids_by_gi(identifiers, key):
    sys.stdout.write("MY VALUES FOR IDENTIFIER {} ARE:\n".format(key))

    return map(lambda x: x[key], identifiers) ## Useful to understand it http://stackoverflow.com/questions/3070242/reduce-python-list-of-objects-to-dict-object-id-object

# http://biopython.org/DIST/docs/api/Bio.PDB.PDBList'-pysrc.html PARA CUALQUIER DUDA EN LOS OBJETOS CONSULTAR
def pdb_download(pdb_id, path=None): # path --> directory
    """ Retrieves a PDB structure file from the PDB server and  stores it in a local file tree. The PDB structure's file name is returned as a single string. If obsolete ``==`` True, the file will be saved in a special file tree.Downloads the structure of the pdb on a file. code = pdb code. """
    sys.stdout.write("DOWNLOADING PDB {}.\n".format(pdb_id))

    pdb_list = PDBList(obsolete_pdb=os.getcwd())
    if path is None:
        file = pdb_list.retrieve_pdb_file(pdb_id)
    else:
        file = pdb_list.retrieve_pdb_file(pdb_id, pdir=path)

    return file

# http://biopython.org/DIST/docs/api/Bio.Align.MultipleSeqAlignment-class.html
def which_msa(msa, input_f, output_f, o_format=None):
    """Runs the MSA with the selected third-party software. TAKE CARE: Muscle only produce FASTA outputs !!!! """
    sys.stdout.write("MSA with {} FROM seq {} to {}\n".format(msa, input_f, output_f))

# Verifying that you have the program installed
    if not find_executable(msa):
        raise UnboundLocalError("THE PROGRAM {} DOES NOT FOUND IN YOUR COMPUTER. ARE YOU SURE IT IS INSTALLED?".format(msa))

# Selecting way of doing the MSA
    if msa == "clustalw":
        command = ClustalwCommandline(msa, infile = input_f, output = o_format, OUTFILE=output_f, type="PROTEIN")
    elif msa == "clustalo":
        command = ClustalOmegaCommandline(msa, infile=input_f, outfile=output_f, verbose=True, auto=True)
    elif msa == "muscle":
        command = MuscleCommandline(msa, input=input_f, out=output_f)
    elif msa == "t_coffee":
        command = TCoffeeCommandline(msa, infile=input_f, output=o_format, outfile=output_f)
    try:
        stdout, stderr = command()

    except ApplicationError:
        raise IOError("THIS INPUT DOES NOT CONTAIN AN ALIGNMENT.")
    return(stdout)


if __name__ == "__main__":
    print("You 're trying to run a program that isn't standalone. This package is part of CMDC.py\n")
