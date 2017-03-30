"""DOCUMENTATION FOR CMDC.py: "CORRELATION MUTATIONS AND DISTANCE CORRELATIONS TO PREDICT AMINOACIDS INTERACTION". \n Written in Python3 by Marcos Camara Donoso & Fernando Pozo Ocampo. Project for the subjects Introduction to Python and Structural Bioinformatics. M.Sc. in Bionformatics for Health Sciences. Pompeu Fabra University.\n\n
distances.py standalone module.\n
Standard modules sys, argparse, os, numpy and BioPython are required in order to execute in a proper way the program.\n
Creating a contact map for a pdb structure given.\n\n
USAGE: $ python3 distances.py pdb_file -atom CA -CA 8 -CB 8 -min 4 """

import sys, argparse, os

import numpy as np

import matplotlib.pyplot as plt

from Bio.PDB.PDBParser import PDBParser

#Adaptation from: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def residue_dist(residue1, residue2, atom):
    """Calculates distances between two residue using the atom selected."""

    sys.stdout.write("Computing calculation of distances between {}-{} and {}-{} at {} atom.\n".format(residue1.get_resname(), residue1.id[1], residue2.get_resname(),residue2.id[1], atom))

    if atom is not None:
        try:
            distance = residue1[atom] - residue2[atom]
        except KeyError:
            sys.stdout.write("Taking C-Alpha for Gly.\n") #Taking into account CA for Gly because it hasn't got CB
            distance = residue1["CA"] - residue2["CA"]
        return distance

#Adaptation from: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def min_dist(residue1, residue2):
    """Calculates the minimum distance between two residues."""

    sys.stdout.write("Computing calculation of minimum distances between {}-{} and {}-{}.\n".format(residue1.get_resname(),residue1.id[1], residue2.get_resname(),residue2.id[1]))

    distances = []

    #Fills the list with all the distances between atoms in two residues and returns the minimum ones.
    for atm1 in residue1:
        for atm2 in residue2:
            distances.append(atm1-atm2)
    return min(distances)

#Adaptation from: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def matrix_dist(chain, atom=None):
    """Calculates the matrix of distances of residues from the same chain."""

    sys.stdout.write("Computing distance matrix for {}\n".format(chain))

    length = len(chain)
    matrix = np.zeros((length, length), np.float)
    #Returns a list of tuples with positions and residues: https://docs.python.org/3/library/functions.html#enumerate
    for row, residue1 in enumerate(chain):
        for col, residue2 in enumerate(chain):
            #Stores the minimun distances or the real distances between two residues taking into account the atom.
            if atom == "min" or atom is None:
                matrix[row, col] = min_dist(residue1, residue2)
            else:
                matrix[row, col] = residue_dist(residue1, residue2, atom)
    sys.stdout.write("Min distance = {}\n".format(np.min(matrix)))
    sys.stdout.write("Max distance = {}\n".format(np.max(matrix)))
    return matrix

def ns_residue_excluder(structure):
    """Exclude non-standard aminoacids of the structure in order to get rid of them"""
    #As we are only interested in standard aminoacids, every kind of a different one will be filtered.
    sys.stdout.write("Parsing structure {}".format(structure))
    residues = tuple(structure.get_residues())
    sys.stdout.write("Total length of structure = {}".format(len(residues)))
    sys.stdout.write("Excluding non-standard aminoacids")
    #Ignores the non-standard and retrieves a new tuple with all the standard residues.
    residues = tuple(filter(lambda x: x.id[0] == " ", residues))
    sys.stdout.write("New length after excluding = {}".format(len(residues)))

    return residues

def contact_map(d_map, atom, dist=None):
    """For a given distance map it selects residues that have significant contact under a distance threshold depending of the atom selected."""

    sys.stdout.write("Filtering atom contacts by taking {}".format(atom))

    if dist is not None:
        distances = dist
    else:
        #Distances retrieved from Adhikari and cheng, 2016.
        distances = dict(CA = 8, CB = 8, min = 4)
        length = len(d_map)
        cmap = np.zeros((length, length), dtype=bool)
        contact = 0
        for i in range(length):
            for k in range(length):
                #An aminoacid can't do contact with itself.
                if abs(i-k) <= 2:
                    pass
                #Threshold of distance for a contact between two aminoacids using the selected atom.
                elif d_map[i][k] <= distances[atom]:
                    contact += 1
                    cmap[i][k] = True
        sys.stdout.write("There are {} contacts".format(contact/2))
        return cmap

#PLOTS FOR CONTACT MAPS
def plot_heatmap(matrix, name_file, title, option, label):
    """Plots a heatmap of a distances."""
    sys.stdout.write("PLOTTING HEATMAP OF DISTANCES... {}".format(name_file))

    #Plots the matrix of distances as a heatmap with its bar.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(title)
    heatmap = plt.pcolormesh(matrix, cmap="plasma")
    legend = plt.colorbar(heatmap)
    legend.set_label(label)
    imgplot = ax.imshow(matrix, interpolation='none')
    #Defines the directory where the plots will be at.
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'plots/')
    #If it doesn't exist it creates the folder.
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    #saves the plot.
    name_out = 'figure5_heatmap_{}_{}.png'.format(name_file, option)
    plt.savefig(results_dir + name_out, format="png")
    return name_out

def plot_contactmap(matrix, name_file, title, option):
    """Plotting classical contact map of the residue"""
    sys.stdout.write("PLOTTING CONTACT MAP FOR {}".format(name_file))

    #Plots the matrix of contacts. Black for contact. White for not.
    fig_b = plt.figure()
    ax = fig_b.add_subplot(111)
    fig_b.suptitle(title)
    imgplot = ax.imshow(matrix, cmap='Greys', interpolation='none')

    #Creates the directory in which the program will save the plots.
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'plots/')

    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    #saves the plot
    name_out = 'figure4_contactmap_{}_{}.png'.format(name_file, option)
    plt.savefig(results_dir + name_out, format="png")
    return name_out


if __name__ == "__main__":

    args = argparse.ArgumentParser(description='Calculates distance map and generates plots', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("file" , help="PDB file.")
    argparser.add_argument("-atom", help = "Residues atom to calculate distances: Alpha carbon, beta carbon and min (minimum distance between atom pairs from each aminoacid).", default="min", choices=["CA", "CB", "min"])
    argparser.add_argument("-CA", help="Threshold for alpha carbon atoms.", type = int, default = 8)
    argparser.add_argument("-CB", help = "Threshold for beta carbons atoms.", type = int, default = 8)
    argparser.add_argument("-min", help = "Set the minimal threshold distance between atoms.", type = int, default = 4)
    arguments = argparser.parse_args()

    sys.stdout = open('stdout_dis.txt', 'w')
    sys.stderr = open('stderr_dis.txt', 'w')

    sys.stdout.write("SESSION HAS ALREADY STARTED WITH THE NEXT PARAMETERS CHOSEN:{}".format(args))
    sys.stdout.write("Analysis of {} using {} atom has begun".format(filename, atom))

    distances = dict( CA = args.CA, CB = args.CB, min = args.min)
    base = os.path.basename(args.file)
    fname = os.path.splitext(base)[0]
    parser = PDBParser(PERMISSIVE=1)

    structure = parser.get_structure("test", args.file)
    residues = ns_residue_excluder(structure)
    dmatrix = dist_matrix(residues, args.atom)
    dtitle = 'HEATMAP OF %s' %(fname)
    #Plotting
    heatmap = plot_heatmap(dmatrix, fname, dtitle, args.atom)
    sys.stdout.write("Heatmap for {} generated ".format(fname))

    cmatrix = cmap(dmatrix, args.atom, distances)
    bintitle = "CONTACT MAP OF %s" %fname
    binmatrix = plot_contactmap(cmatrix, fname, bintitle, atom)

    sys.stdout.write("Contact map plot for %s generated".format(fname))
    sys.stdout.write("Contact maps computed. Exiting program...")

    exit(0)
    sys.stdout.close()
    sys.stderr.close()
