#!/usr/bin/python3
"""Blah
"""

__version__ = "0.0.0"

import argparse
import pandas as pd
from straw import straw
from joblib import Parallel, delayed
#import juice


def get_contacts_frame(optArgs, chrA, chrB):
    """Extract the confrom numba import njit, prangetact matrix in flat format from the .hic file and return a double indexed data frame with the contacts.

    """
    chrA = "chr" + chrA
    chrB = "chr" + chrB
    res = straw(optArgs.norm, optArgs.infile, chrA, chrB, optArgs.type, optArgs.binSize)
    multi_index = pd.MultiIndex.from_tuples(tuples = list(zip(res[0], res[1])), names = ['start', 'stop'])
    cont = "contacts_" + chrA + "*" + chrB
    dc = pd.DataFrame(data = {cont:res[2]}, index = multi_index)
    return(dc)


def extract_contacts(optArgs, chrA, chromosomes):
    """The paralelisation wrapper of the whole extract process.

    Works one chromosome at a time (i.e. it paralelises all the contacts beteen a given (chrA) and all the rest of the chromosoems.)
    """
    all_chrom_res = Parallel(n_jobs=-1, backend = "multiprocessing")(delayed(get_contacts_frame)(optArgs, chrA, chrB) for chrB in chromosomes)
    return(dict(zip(chromosomes, all_chrom_res)))


parser = argparse.ArgumentParser(prog='pyna_collada', description='Blah blah', epilog="Authors: Costas Bouyioukos, 2019-2020, Universite de Paris et UMR7216.")
parser.add_argument('infile', type=str, metavar="input_file", help='Filename (or path) of a hic file (NO option for STDIN).')
#parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path to output FASTA file. (or STDOUT).")
parser.add_argument('-n', '--normalisation', nargs="?", default='VC', metavar="Norm. meth.", type=str, help="Choise of a normalisation method from the Juice suite or straw (Default: VC).", dest="norm")
parser.add_argument('-t', '--type', nargs="?", default='BP', metavar="Type", type=str, help="Choise of a measure (Default: BP).", dest="type")
parser.add_argument('-b', '--bin-size', help="Seelction of the bin size of the hi-c map (i.e. resolution). (Default=25000).", type=int, default=25000, dest="binSize", metavar="Bin Size")
parser.add_argument(-g, --gene-list, type=argparse.FileType('r'), default=None, dest="genesCoord", metavar="gene's coords", help="A list of genes (or genomic locations) of interest and their genomic coordinates. The full length of gene is considered here.")
parser.add_argument('-v', '--version', action='version', version='%(prog)s  v. {version}'.format(version=__version__))
#TODO fix the argument ranges of accepted values form straw.
#TODO arguments: Add argument for chromosomes/organism.
#

# Parse the command line arguments.
optArgs = parser.parse_args()

chromosomes=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

chromosomes=["1","2"]

# Read the file and put results in a dictionary -> data frame structure.
contacts = {}
for chrA in chromosomes:
    contacts[chrA] = extract_contacts(optArgs, chrA, chromosomes)
