#!/usr/bin/python3
"""Blah
"""

__version__ = "0.0.0"

import argparse
import pandas as pd
from straw import straw
#import parallel
#import juice


parser = argparse.ArgumentParser(prog='pyna_collada', description='Blah blah', epilog="Authors: Costas Bouyioukos, 2019-2020, Universite de Paris et UMR7216.")
parser.add_argument('infile', type=str, metavar="input_file", help='Filename (or path) of a hic file (NO option for STDIN).')
#parser.add_argument("outfile", nargs='?', default='-', type=argparse.FileType('w'), metavar='output_file', help="Path to output FASTA file. (or STDOUT).")
parser.add_argument('-n', '--normalisation', nargs="?", default='VC', metavar="Norm. meth.", type=str, help="Choise of a normalisation method from the Juice suite or straw (Default: VC).", dest="norm")
parser.add_argument('-t', '--type', nargs="?", default='BP', metavar="Type", type=str, help="Choise of a measure (Default: BP).", dest="type")
parser.add_argument('-b', '--bin-size', help="Seelction of the bin size of the hi-c map (i.e. resolution). (Default=25000).", type=int, default=25000, dest="binSize", metavar="Bin Size")
parser.add_argument('-v', '--version', action='version', version='%(prog)s  v. {version}'.format(version=__version__))
#TODO fix the argument ranges of accepted values form straw.

# Parse the command line arguments.
optArgs = parser.parse_args()

# Read the file and put results in a pandas data frame.
res = straw(optArgs.norm, optArgs.infile, "chr1", "chr2", optArgs.type, optArgs.binSize)
dd = pd.DataFrame(data = {"START":res[0], "STOP":res[1], "CONTACT":res[2]})
print(dd.head())
