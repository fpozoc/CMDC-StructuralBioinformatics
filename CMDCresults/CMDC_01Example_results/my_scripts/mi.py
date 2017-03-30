""""DOCUMENTATION FOR CMDC.py: "CORRELATION MUTATIONS AND DISTANCE CORRELATIONS TO PREDICT AMINOACIDS INTERACTION".\n
Written in Python3 by Marcos Camara Donoso & Fernando Pozo Ocampo. Project for the subjects Introduction to Python and Structural Bioinformatics. M.Sc. in Bionformatics for Health Sciences. Pompeu Fabra University.\n\n
mi.py standalone module.\n
Standard modules sys, argparse, os, copy, math, numpy, pylab, pytplot and BioPython are required in order to execute in a proper way the program.\n
This script calculates mutual information normalized Z-scores between all positions in a single protein, make a correction removing possible noise in the study of MSA and returns the correlations for all possible pair of residues.


\n\n\n USAGE: $ python3 mi.py msa_file.fa target_pdbcode -threshold 2.0 -b 20 -gaps -low 0.3 -high 0.9"""

import sys, argparse, copy, math, os

import matplotlib.pyplot as plt

import numpy as np

from Bio import AlignIO

from Bio.PDB import PDBList

import pylab as plot
import matplotlib.pyplot as plt
import numpy, scipy, pylab
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator

## This paper describes the next functions about calculous of Mutual Information Bremm, S., Schreck, T., Boba, P., Held, S., Hamacher, K.: Computing and visually analyzing mutual information in molecular co-evolution.
def frequencies_msa(msa_out, column):
    """Calculate frequencies of each column of multiple sequence alignment."""
    sys.stdout.write("FRECUENCIES OF MSA !\n")

    # Complete set of aminacids
    aminoacids = set('ACDEFGHIKLMNPQRSTVWY-')
    frequencies = dict.fromkeys(aminoacids, 0)

    # Necesitas el índice de cada columna de aminoácidos !
    column_msa = msa_out[:, column] # alineamiento que le das de entrada e índices
    for letter in column_msa:
        try:
            frequencies[letter] += 1
        except KeyError:
            pass
    for k in frequencies:
        frequencies[k] = frequencies[k]/len(column_msa)

    return frequencies

def jfrequencies_msa(msa_out, column1, column2):
    """Calculate frequencies of pair of column of multiple sequence alignment."""
    sys.stdout.write("JOINT FREQUENCY OF {} & {}\n".format(column1, column2))

    # Vamos a hacer lo mismo que en la función anterior pero esta vez guardando un tuple de residuos para cada par de columnas
    aminoacids = set('ACDEFGHIKLMNPQRSTVWY-')

    # Creating a set of indexes
    index_set = set()
    for i in aminoacids:
        for j in aminoacids:
            # Adding a tuple to the set!
            index_set.add((i, j))

    # Cuidado con el uso correcto de fromkeys para crear el nuevo diccionario https://www.tutorialspoint.com/python/dictionary_fromkeys.htm
    frequencies = dict.fromkeys(index_set, 0)
    column_x = msa_out[:, column1]
    column_y = msa_out[:, column2]

    # Para pasar por todas las posiciones del msa_out
    for i in range(len(msa_out)):
        try:
            frequencies[(column_x[i], column_y[i])] += 1
        except KeyError:
            pass
    # Same for loop
    for k in frequencies:
        frequencies[k] = frequencies[k]/len(msa_out)

    return frequencies # Returning join frequency

# Calculating entropy
def entropy_msa(msa_out, column, logarithm_base):
    """Calculate ENTROPIES of each column of multiple sequence alignment !!! """
    sys.stdout.write("ENTROPY OF {} COLUMN\n".format(column))

    frequencies = frequencies_msa(msa_out, column)
    entropy_msa = 0
    for k in frequencies:
        if frequencies[k] != 0:
            entropy_msa -= frequencies[k]*math.log(frequencies[k], logarithm_base)
    return entropy_msa

# The same as before with the calculous of entropy (joint) this time
def joint_entropy_msa(msa_out, column1, column2, logarithm_base):
    """Calculate JOINT ENTROPY ENTROPIES of each column of multiple sequence alignment !!! """
    sys.stdout.write("JOINT ENTROPY OF {} & {}\n".format(column1, column2))

    joint_frequencies = jfrequencies_msa(msa_out, column1, column2)
    jentropy_msa = 0
    for k in joint_frequencies:
        if joint_frequencies[k] != 0:
            jentropy_msa -= joint_frequencies[k]*math.log(joint_frequencies[k], logarithm_base)

    return jentropy_msa

def minmax_entropy(msa_out, entropy_min, entropy_max, logarithm_base):
    """Minimum and maximum columns entropies."""

    sys.stdout.write("MIN & MAX COLUMNS ENTROPY ABOVE THRESHOLD.\n")

    minlist = []
    maxlist = []

    #http://biopython.org/DIST/docs/api/Bio.Align.Generic-pysrc.html#Alignment.get_alignment_length
    for i in range(msa_out.get_alignment_length()):
        if entropy_msa(msa_out, i, logarithm_base) <= entropy_min:
            minlist.append(i)
        if entropy_msa(msa_out, i, logarithm_base) >= entropy_max:
            maxlist.append(i)

    return (minlist, maxlist)


def mutual_information(msa_out, column1, column2, logarithm_base=20): # Base of the logarithm
    """Calculates Mututal Information with the formula"""
    sys.stdout.write("MUTUAL INFminmax_entropyn".format(column1, column2))

    entropy_1 = entropy_msa(msa_out, column1, logarithm_base)
    entropy_2 = entropy_msa(msa_out, column2, logarithm_base)

    return entropy_1 + entropy_2 - joint_entropy_msa(msa_out, column1, column2, logarithm_base)


def mutual_information_matrix(msa_out, logarithm_base=20):
    """Calculates an Square Mutual Information Matrix"""
    sys.stdout.write("MUTUAL INFORMATION MATRIX\n")

    seqs = msa_out.get_alignment_length()

    #https://docs.scipy.org/doc/numpy/reference/generated/numpy.empty.html
    matrix = np.empty([seqs, seqs], dtype=float, order='F') # Takes empty matrix

    # Filling the empty matrix created
    for j in range(seqs):
        for i in range(seqs):
            if j <= i:
                matrix[i, j] = mutual_information(msa_out, i, j, logarithm_base)
            else:
                matrix[i, j] = copy.copy(matrix[j, i])

    return matrix


def matrix_hits(binary_matrix):
    """Gets values 1 in binary matrix."""
    sys.stdout.write("CONTACTS PREDICTED BELOW...\n")

    counter = 0
    (n1, n2) = binary_matrix.shape

    for i in range(n1):
        for j in range(n2):
            if j > i:
                counter += binary_matrix[i, j]

    return counter

# Computing with numpy mean and standard deviation of our values and creating Z-scores
def standardise_matrix(mz):
    """Normalizes the matrix with Z-Scores."""
    sys.stdout.write("MUTUAL INFORMATION VALUES TO Z-SCORES...\n")

    array = []
    (n1, n2) = mz.shape
    for i in range(n1):
        for j in range(n2):
            if i != j:
                array.append(mz[i, j])

    #  https://docs.scipy.org/doc/numpy/reference/generated/numpy.mean.html
    mean = np.mean(array)
    std_dev = np.std(array)
    matrix = np.empty([n1, n2], dtype=float, order='F')

    for i in range(n1):
        for j in range(n2):
            if i == j:
                matrix[i, j] = 0
            else:
                matrix[i, j] = (mz[i, j]-mean)/std_dev

    return matrix


####### Corrección implementada extraída de la publicación en 2009 de Byung-Chul y Dongsup en la que confirma que esta corrección elimina ruido de fondo  en el cálculo de correlaciones en la coevolución de residuo
####### "A new method for revealing correlated mutations under thestructural and functional constraints in proteins"
####### Primero se calcula fórmula para CPS y luego se normaliza el valor (NCPS)
def CPS(msa_out, column1, column2, logarithm_base=20):
    """Calculates coevolutionary pattern similarity."""
    sys.stdout.write("CO-EVOLUTIONARY PATTERN SIMILARITY(CPS)\n")

    seqs = msa_out.get_alignment_length()

    if seqs <= 1:
        raise ValueError("You need more sequences in your multiple sequence alignment.\n")

    mim = mutual_information_matrix(msa_out, logarithm_base)
    cps = 0

    ## Aquí tenemos la primera fórmula
    for k in range(seqs):
        if k != column1 and k != column2: # if k ne i,j valores de la diagonal de la columna
            cps += mim[column1, k]*mim[k, column2]
    cps = cps/(seqs-2)

    return cps

###### Lo mismo pero para normalizar el valor anterior
def NCPS_matrix(msa_out, logarithm_base=20):
    """Calculates normalized coevolutionary pattern similarity"""
    sys.stdout.write("NORMALIZED CO-EVOLUTIONARY PATTERN SIMILARITY(NCPS)\n")

    seqs = msa_out.get_alignment_length()

    if seqs <= 1:
        raise ValueError("You need more sequences in your multiple sequence alignment. Lenght -> {}\n".format(seqs))

    mim = mutual_information_matrix(msa_out, logarithm_base)
    cps_matrix = np.empty([seqs, seqs], dtype=float, order='F')
    all_m = 0

    for j in range(seqs):
        for i in range(seqs):
            cps = 0
            if j <= i:
                for k in range(seqs):
                    if k != i and k != j:
                        cps += mim[i, k]*mim[k, j]
                cps_matrix[i, j] = cps / (seqs-2)
            else:
                cps_matrix[i, j] = copy.copy(cps_matrix[j, i])
            all_m += cps_matrix[i, j]
    all_m = all_m/(seqs*(seqs-1))
    all_m = math.sqrt(all_m)

    return cps_matrix * 1/all_m

## Vamos a sacar el MIc value que calcula en el paper de NCPS
def MIc(msa_out):
    """Gets MI-NCPS"""
    sys.stdout.write("CALCULATES MI-NCPS AND NORMALIZED VALUES...\n")

    MIc = mutual_information_matrix(msa_out) - NCPS_matrix(msa_out)

    return (MIc, standardise_matrix(MIc))


def prior_position(position, deleted_position):
    """Prior positions in Multiple Sequence Alignment."""
    sys.stdout.write("ORIGINAL POSITIONS IN THE MATRIX.\n")

    difference = 0
    deleted_pos = sorted(deleted_position)
    for n in range(len(deleted_position)):
        if deleted_position[n] <= position + difference:
            difference += 1
        else:
            break

    return position + difference

# IMPORTANT !
def binary_matrix(matrix, level):
    """Binary matrix with positions above threshold."""
    sys.stdout.write("PREDICTING RESIDUES PAIR CONTACTS.\n")

    (n1, n2) = matrix.shape
    out_matrix = np.empty([n1, n2], dtype=float, order='F')

    for i in range(n1):
        for j in range(n2):
            if i == j:
                out_matrix[i, j] = 0
            elif matrix[i, j] >= level:
                out_matrix[i, j] = 1
            else:
                out_matrix[i, j] = 0

    return out_matrix

def prior_index_position(binary_matrix, gap_list, extreme_list):
    """Prior positions indexes in Multiple Sequence Alignment. Retrieves that indexes with value 1 in previous binary matrix."""
    sys.stdout.write("ORIGINAL COORDINATES OF THE MATRIX!\n")

    pair_residues = set()
    (n1, n2) = binary_matrix.shape
    for i in range(n1):
        i_reconstruct = prior_position(i, extreme_list)
        for j in range(n2):
            if j > i and binary_matrix[i, j] == 1:
                j_reconstruct = prior_position(j, extreme_list)
                residue1 = prior_position(i_reconstruct, gap_list)
                residue2 = prior_position(j_reconstruct, gap_list)
                pair_residues.add(tuple(sorted([residue1, residue2])))

    return pair_residues


def original_positions(matrix, gap_list, extreme_list):
    """Retrieves a dictionary with keys given by the original coordinates
    of the set of pairs of residues corresponding to all spots in matrix
    above the diagonal, and values corresponding to the values in matrix."""

    sys.stdout.write("Retrieving all positions with original coordinates.\n")

    pair_residues_dict = {}
    (n1, n2) = matrix.shape
    for i in range(n1):
        i_reconstruct = prior_position(i, extreme_list)
        for j in range(n2):
            if j > i:
                j_reconstruct = prior_position(j, extreme_list)
                residue1 = prior_position(i_reconstruct, gap_list)
                residue2 = prior_position(j_reconstruct, gap_list)
                value = matrix[i, j]
                pair_residues_dict[tuple(sorted([residue1, residue2]))] = value
    return pair_residues_dict

def distances_matrix_MIvalues(pair_residues_dict, dist_matrix):
    """Distances matrix compare Mutual Information Values."""
    sys.stdout.write("DISTANCES MATRIX - MI VALUES.\n")

    # IMPORTANTE PRINT DE CONTROL
    print(pair_residues_dict)
    list_coor = pair_residues_dict.keys()
    list_mutual_information = list(pair_residues_dict.values())
    list_dist = []

    for element in list_coor:
        list_dist.append(dist_matrix[element[0]][element[1]])
    if len(list_dist) == len(list_mutual_information):

        return (list_dist, list_mutual_information)
    else:
        raise ValueError("PROBLEM WITH LISTS!\n")

##### PRUNING COLUMNS. PLEASE, CHECK REFERENCE
##### Reliable and robust detection of coevolving protein residues. Chan Seok 2012
def remove_columns(msa_out, columns_list):
    """Removing corresponding columns from Multiple sequence alignment."""
    sys.stdout.write("PRUNING...\n")

    columns_list = sorted(columns_list)

    for i in range(len(columns_list)):
        msa_out = msa_out[:, : columns_list[i]-i] + msa_out[:, columns_list[i]-i+1:]
    return msa_out


def remove_columns_gaps(msa_out, identifier):
    """Removing gaps columns from Multiple sequence alignment."""

    sys.stdout.write("PRUNING GAPS...\n")

    for i in range(len(msa_out)):
        if msa_out[i].id == identifier:
            row = i
    columns_list = []

    for j in range(msa_out.get_alignment_length()):
        if msa_out[row].seq[j] == '-':
            columns_list.append(j)

    return remove_columns(msa_out, columns_list)

def get_columns_onegap(msa_out):
    """Gets one gaps columns."""
    sys.stdout.write("LOOGKING FOR POSITIONS WITH GAPS.\n")

    columns_list = []

    for j in range(msa_out.get_alignment_length()):
        for i in range(len(msa_out)):
            if msa_out[i].seq[j] == '-':
                columns_list.append(j)
                break

    return columns_list


######## PLOTS ########
### Revisar Correlated Mutations and Residue Contacts in Proteins. Valencia 1994

def plot_MI_res(zmatrix, title, name_file, option):
    """Calculates parameters and plots MI z-score for predicted contacts per column vs residues"""
    sys.stdout.write("PLOTTING MI Z-SCORE FOR PREDICTED CONTACTS OF {}\n".format(name_file))

    array_max = np.amax(zmatrix, axis=1)
    array_mean = np.mean(zmatrix, axis=1)
    array_sd = np.std(zmatrix, axis=1)
    array_sd_mean = array_mean + array_sd
    pos = range(len(array_max))

    fig = plt.figure()
    fig.suptitle(title)
    ax = plt.subplot(111)

    ax.xaxis.set_major_locator(MaxNLocator(10))
    ax.yaxis.set_major_locator(MaxNLocator(7))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))

    plt.plot(pos, array_max, 'og-', label='Max Z-score MI',  linewidth=1)
    plt.plot(pos, array_mean, 'ob-', label='Mean Z-score MI',  linewidth=1)
    plt.plot(pos, array_sd_mean,'r-', label='SD area', linewidth=1)
    ax.fill_between(pos, array_mean, array_sd_mean, color='red', alpha=0.5)
    plt.legend( loc='upper right', fancybox=True, framealpha=0, borderaxespad=0., prop={'size':8})

    plt.xlabel('Residues', fontsize=10)
    plt.ylabel('Z-Score MI', fontsize=10)

    for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(0.5)

    plt.subplots_adjust(bottom=0.18)
    plt.subplots_adjust(left=0.14)
    plt.subplots_adjust(right=0.96)

    plt.tick_params(axis='both', which='major', labelsize=8)
    plt.tick_params(axis='both', which='minor', labelsize=0)

    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'plots/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    name_out = "figure1_Z-scoreMI_max_mean_{}_{}".format(name_file, option)

    plt.savefig(results_dir + name_out ,format="png")
    return name_out

def plot_superimp_matrix(matrix, matrix2, name_file, title, label, option):
    """Predicted contacts with Z-score exceeding 2 over the rest of the Mutual Information Values"""
    sys.stdout.write("PLOTTING SUPERIMPOSE MATRIX OF {}\n".format(name_file))

    fig_b = plt.figure()
    ax = fig_b.add_subplot(111)
    fig_b.suptitle(title)
    heatmap = plt.pcolormesh(matrix2, cmap = 'inferno', alpha=0.3)
    legend = plt.colorbar(heatmap)
    legend.set_label(label)
    ax.imshow(matrix2, interpolation='none')
    ax.imshow(matrix, cmap='Greys', interpolation='none')

    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'plots/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    name_out = 'figure2_supermatrix_{}_{}.png'.format(name_file, option)
    plt.savefig(results_dir + name_out, format="png")
    return name_out

def plot_accuracy(cutoff_list, hit_list, precision_list, name_file):
    """Accuracy, predicted contacts and correlation cutoff rc"""
    sys.stdout.write("PLOTTING ACCURACY...\n")

    fig, ax1 = plt.subplots()
    ax1.plot(cutoff_list, hit_list, 'og-')
    ax1.set_xlabel('correlation cutoff rc')
    ax1.set_ylabel('number of predicted contacts', color='g')
    ax2 = ax1.twinx()
    ax2.plot(cutoff_list, precision_list, '^r-')
    ax2.set_ylabel('accuracy', color='r')

    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'plots/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    name_out = 'figure6_accuracy_cutoffs_{}.png'.format(name_file)
    plt.savefig(results_dir + name_out, format="png")

    return name_out

def plot_scatter_distance(list_dist, list_mutual_information, name_file):
    """Scatter plot of residue pair distances"""
    sys.stdout.write("SCATTER PLOTTING...\n")

    fig = plt.figure()
    fig_title = 'Z-score MI and DISTANCE BETWEEN RESIDUE PAIRS OF {}'
    fig.suptitle(fig_title.format(name_file))
    assert len(list_dist) == len(list_mutual_information)
    plt.scatter(list_mutual_information, list_dist, color='blue', s=5, edgecolor='none')
    plt.xlabel('correlation (Z-scoresMI)', fontsize=16)
    plt.ylabel('Distance (Angstroms)', fontsize=16)

    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'plots/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    file_out = 'figure3_scatterplot_distances_{}.png'.format(name_file)
    plt.savefig(results_dir + file_out, format="png")
    return file_out

def precision_analysis(mutual_information_val, cont_m, gapped_list, mm, l=0.0, h=3.0, num=60):
    """Calculating accuracy (Precision of analysis) and hits for diffent threshhold levels of MI values."""
    sys.stdout.write("CALCULATING ACCURACY...\n")

    cutoff_list = []
    hit_list = []
    precision_list = []

    for cutoff in np.linspace(l, h, num):
        cutoff_list.append(cutoff)
        tmatrix = binary_matrix(mutual_information_val, cutoff)
        hits = matrix_hits(tmatrix)
        hit_list.append(hits)
        distances_residue_pairs = prior_index_position(tmatrix, gapped_list, mm)

        count = 0
        for rp in distances_residue_pairs:
            count += cont_m[rp[0], rp[1]]
        precision_list.append(count/hits)
    else:

        return(cutoff_list, hit_list, precision_list)

if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description="Computations of mutual information", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # compulsory input
    argparser.add_argument("i", help="Multiple Sequence Alignment FASTA file")
    argparser.add_argument("id", help="Target identifier")
    # optional input
    argparser.add_argument("-output",
                           help="""Residues predicted output""",
                           default="my_predicted_residuecontacts.out")
    argparser.add_argument("-threshold",
                           help="""Threshold level to make predictions.""",
                           type=float,
                           default=2.0)
    argparser.add_argument("-b",
                           help="""Base of the logarithms (entropy and mutual information in mi.py are used for that)""",
                           type=int,
                           default=20)
    argparser.add_argument("-gaps", help="""If user selects, then remove those columns of the MSA which have at least one gap.""",
                           action='store_true',
                           default=False)
    argparser.add_argument("-low",
                           help="""Mininum entropy threshold allowed for each column in the MSA.""",
                           type=float,
                           default=0.3)
    argparser.add_argument("-high",
                           help="""Maximum entropy threshold allowed for each column in the MSA.""",
                           type=float,
                           default=0.9)
    args = argparser.parse_args()


    sys.stdout = open('stdout_MI.txt', 'w')
    sys.stderr = open('stderr_MI.txt', 'w')

    sys.stdout.write("SESSION HAS ALREADY STARTED WITH THE NEXT PARAMETERS CHOSEN:{}".format(args))

    # IMPORTANTE QUE EL MSA ESTÉ EN FASTA
    alignment = AlignIO.read(args.i, "fasta")

    edited = remove_columns_gaps(alignment, args.id)
    gapped_list = []

    if args.gaps:
        gapped_list = get_columns_onegap(edited)
    edited = remove_columns(edited, gapped_list)
    (minlist, maxlist) = minmax_entropy(edited, args.low, args.high, args.b)
    edited = remove_columns(edited, minlist+maxlist)

    MI_matrix = mutual_information_matrix(edited, args.b)
    ncps_array = NCPS_matrix(edited, args.b)
    MIc_matrix = MI_matrix - ncps_array
    mutual_information_matrix = standardise_matrix(MIc_matrix)

    title_mutualinfo_scores = 'MI INFORMATION VALUES OF {}'.format(args.i)
    plot_heatmap(mutual_information_matrix, args.i, title_mutualinfo_scores, args.low, "zMIc")
    tmatrix = binary_matrix(mutual_information_matrix, args.threshold)
    title_mutualinfo_scores_b = "MI VALUES PREDICTED CONTACTS WITH Z-score MI > {} OF {}".format(args.threshold, args.i)
    plot_superimp_matrix(tmatrix, args.i, title_mutualinfo_scores_b, args.threshold)
    distances_residue_pairs = prior_index_position(tmatrix, gapped_list, minlist + maxlist)
    fd = open(args.output, "w")
    for residue_pair in sorted(distances_residue_pairs):
        fd.write("%d %d\n" % residue_pair)
    fd.close()

    sys.stdout.write("Finally...\n")
    exit(0)

    sys.stdout.close()
    sys.stderr.close()
