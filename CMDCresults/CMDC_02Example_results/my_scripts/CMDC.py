"""DOCUMENTATION FOR CMDC.py: "CORRELATION MUTATIONS AND DISTANCE CORRELATIONS TO PREDICT AMINOACIDS INTERACTION". \nWritten in Python3 by Marcos Camara Donoso & Fernando Pozo Ocampo. Project for the subjects Introduction to Python and Structural Bioinformatics. M.Sc. in Bionformatics for Health Sciences. Pompeu Fabra University.\n\n
THIS IS THE MAIN SCRIPT OF THE PROGRAM WORKFLOW.\n

Standard modules sys, argparse, os, BioPython are required in order to execute in a proper way the program. Local and Non-standard modules extract_sequences.py, mi.py and distances.py are also neccesary.\n\n

The user must introduce a pdb identifier of a protein. Then, a BLAST alignment of that sequence will be provided. Retrieving each sequence of this alignment from PDB the program will return a MSA alignment with its program selected. Consequently, it will calculates the right correlations between residues with the right information theory calculus. At the same time it will relacionate that scores calculating distances between aminoacids and stablishing possible contacts. Correlations between this possible contacts, distances and correlated mutations will be plotted in order to show graphically to the user the result in its protein."""

import sys, argparse, os

from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio import SeqIO, Entrez, AlignIO
from Bio.Data import SCOPData

import extract_sequences as extseq
import mi
import distances as dis

Entrez.tool = "CMDC_project.py"  # http://biopython.org/DIST/docs/api/Bio.Entrez-module.html
Entrez.email = "fernando.pozo01@estudiant.upf.edu"

# CMDC.py
if __name__ == '__main__':
    default_help = argparse.ArgumentDefaultsHelpFormatter
    inp = argparse.ArgumentParser(description ='Make a structure of a real pdb', formatter_class=default_help)
    inp.add_argument("input", help="Target sequence id or filename")
    # argument options for distances.py module functions
    inp.add_argument("-atom", help="""Residues atom to calculate distances: Alpha carbon, beta carbon and min (minimum distance between atom pairs from each aminoacid).""", choices=["CA", "CB", "min"], default="min")
    inp.add_argument("-CA", help="""Threshold for alpha carbon atoms.""", type=int, default=8)
    inp.add_argument("-CB", help="""Threshold for beta carbons atoms.""", type=int, default=8)
    inp.add_argument("-min", help="""Set the minimal threshold distance between atoms.""", type=int, default=4)

    # argument options for extract_sequences.py module functions
    inp.add_argument("-blast", help="""Type of BLAST.""", default="blastp") # It can be implement another type of BLAST in the future develop
    inp.add_argument("-db", help="""Which database do you prefer?""", choices=["pdb", "swissprot", "nr"], default="nr")
    inp.add_argument("-seqs", help="""BLAST hits selected.""", type=int, default=200)
    inp.add_argument("-filt", help="""If present, then don't filter the BLAST output by genus for attaining non-redundancy; otherwise filter by genus.""", action='store_false', default=True)

    # argument options for mi.py module functions
    inp.add_argument("-gaps", help="""If user selects, then remove those columns of the MSA which have at least one gap.""", action='store_true', default=False)

    # argument options for mi.py module functions
    inp.add_argument("-b", help="""Base of the logarithms (entropy and mutual information in mi.py are used for that)""", type=int, default=20)
    inp.add_argument("-low", help="""Mininum entropy threshold allowed for each column in the MSA.""", type=float, default=0.3)
    inp.add_argument("-high", help="""Maximum entropy threshold allowed for each column in the MSA.""", type=float, default=0.9)
    inp.add_argument("-msa", help="""Multiple Sequence Alignment method""", default="clustalw", choices=["clustalw", "clustalo", "muscle", "t_coffee"])

    args = inp.parse_args()
    sys.stdout = open('stdout.txt', 'w')
    sys.stderr = open('stderr.txt', 'w')

    sys.stdout.write("SESSION HAS ALREADY STARTED WITH THE NEXT PARAMETERS CHOSEN:{}\n\n".format(args))

    try:
        args.input
    except AttributeError:
        print ("Are you sure your input sequence is right?!")
    else:
        # POR SI ACASO AÑADIMOS EL PERMISIVE WAY http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
        parser = PDBParser(PERMISSIVE=1)
        if os.path.isfile(args.input):
            pdbpath = args.input
        else:
            pdbl = PDBList()
            try:
                pdbpath = extseq.pdb_download(args.input, os.getcwd())
            except:
                raise FileExistsError("Is your query format right?!")
        structure = parser.get_structure(args.input+"_pdb_query", pdbpath)
        residues = dis.ns_residue_excluder(structure)
        s = ""
        for residue in residues:
            s += SCOPData.protein_letters_3to1.get(residue.get_resname(), 'X')
        seq = Seq(s, generic_protein)

        # Compute distances and contact between residues
        dist_matrix = dis.matrix_dist(residues, args.atom)
        cont_matrix = dis.contact_map(dist_matrix, args.atom)


        # Set up the files for calling BLAST

        blast_query_name = args.input+"_blast_query.fa"
        blast_out_name = args.input+"_blast.out"

        fd = open(blast_query_name, "w")
        fd.write(">%s\n%s" % (args.input+"_blast_query", seq))
        fd.close()


        # Enter the query as first element in the output file
        file_out = open(blast_out_name, "w")
        record = SeqIO.read(args.input+"_blast_query.fa", format="fasta")
        SeqIO.write(record, file_out, "fasta")

        # Call run_BLAST and write the output in the output file
        # So far it is not working
        blast_result = extseq.BLAST_online(blast_query_name, args.blast, args.db, args.seqs)
        ides = extseq.BLAST_filter(blast_result, args.filt)
        ids = list(extseq.filt_ids_by_gi(ides, "gi"))

        ## IF YOU WANNA TRY TO MAKE THE STUDY MORE RELIABLE IS BETTER TO IMPROVE THE NUMBER OF MINIMUM HITS. AT LEAST 20
        if len(ids) <= 3:
            msg_err = "FEW HITS ON BLAST...SORRY!"
            raise ValueError(msg_err)
        print(len(ids), ids)
        SeqIO.write(extseq.extract_sequence_from_BLASTids(ids), file_out, "fasta")
        file_out.close()

        # See extract_sequences.py for check the variables !
        extseq.which_msa(args.msa, blast_out_name, args.input+".fa", "fasta")
        alignment = AlignIO.read(args.input+".fa", "fasta")


        # CALCULOUS AND COMPUTATIONS OF FUNCTIONS
        edited = mi.remove_columns_gaps(alignment, args.input+"_blast_query")
        gapped_list = []
        if args.gaps:
            gapped_list = mi.get_columns_onegap(edited)
        edited = mi.remove_columns(edited, gapped_list)
        (minlist, maxlist) = mi.minmax_entropy(edited, args.low, args.high, args.b)
        mm = minlist + maxlist
        edited = mi.remove_columns(edited, mm)
        MI_matrix = mi.mutual_information_matrix(edited, args.b)
        ncps_array = mi.NCPS_matrix(edited, args.b)
        MIc_matrix = MI_matrix - ncps_array
        mutual_information_matrix = mi.standardise_matrix(MIc_matrix)

        # PLOTS AND OUTPUTS
        title_dist = 'HEATMAP OF {}'.format(args.input)
        dis.plot_heatmap(dist_matrix, args.input+"_d", title_dist, args.atom, "Angstroms")
        title_binary = 'CONTACT MAP OF {}'.format(args.input)
        dis.plot_contactmap(cont_matrix, args.input+"_c", title_binary, args.atom)
        tmatrix = mi.binary_matrix(mutual_information_matrix, 2)
        std_cont = mi.prior_index_position(tmatrix, gapped_list, mm)
        name_out = "my_predicted_residuecontacts_{}.out"
        with open(name_out.format(args.input), "w") as out_f:
            out_f.write(repr(std_cont))
        title_mutualinfo_scores_b = "PREDICTED CONTACTS WITH MI Z-scores > 2 (black) OVER THE REST OF {}"
        mi.plot_superimp_matrix(tmatrix, mutual_information_matrix, args.input+"_p", title_mutualinfo_scores_b.format(args.input), "Z-score MI", args.atom)
        mutual_information_distances = mi.distances_matrix_MIvalues(mi.original_positions(mutual_information_matrix, gapped_list, mm), dist_matrix)
        mi.plot_scatter_distance(mutual_information_distances[0], mutual_information_distances[1], args.input)
        title_mutualinfo_scores_res = "Z-score MI PER PREDICTED CONTACT RESIDUES OF {}".format(args.input)
        mi.plot_MI_res(mutual_information_matrix, title_mutualinfo_scores_res, args.input+"_p", args.atom)
        cutoff_list, hit_list, precision_list = mi.precision_analysis(mutual_information_matrix, cont_matrix, gapped_list, mm, 0.0, 3.0, 60)
        mi.plot_accuracy(cutoff_list, hit_list, precision_list, args.input)

        sys.stdout.close()
        sys.stderr.close()
