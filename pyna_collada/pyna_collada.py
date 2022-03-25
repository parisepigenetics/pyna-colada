#!/usr/bin/python3
"""Pyna collada is a tool that draws chromosome contacts of specific prespefied regions from a whole genome Hi-C .hic file chromosome map. Visulaises the resulting contact sub-matrices in interective html files.
"""

__version__ = "0.5"

# import sys
import argparse
import math
import os.path
import numpy as np
import pandas as pd
import multiprocessing as mp
import plotly
import plotly.graph_objects as go
from hicstraw import straw


def get_contacts_frame(optArgs, chrA, chrB):
    """Extract the contact matrix in flat format from the .hic file and return a double indexed data frame with the contacts.
    """
    chrAs = "Chr" + chrA
    chrBs = "Chr" + chrB
    #print(f'chr1 {chrAs} chr2 {chrBs}')
    res = straw("observed", optArgs.norm, optArgs.infile, chrA, chrB, optArgs.type, optArgs.binSize)
    multi_index = pd.MultiIndex.from_tuples(tuples=list(zip(res[0], res[1])), names=["start", "stop"])
    cont = "contacts_" + chrA + "*" + chrB
    dc = pd.DataFrame(data={cont: res.counts}, index=multi_index)
    return dc


def extract_contacts(optArgs, chrA, chromosomes):
    """The paralelisation wrapper of the whole extract process.

    Works one chromosome at a time (i.e. it paralelises all the contacts beteen a given (chrA) and all the rest of the chromosomes.)
    """
    all_chrom_res = Parallel(n_jobs=-1, backend="multiprocessing")(delayed(get_contacts_frame)(optArgs, chrA, chrB) for chrB in chromosomes)
    return dict(zip(chromosomes, all_chrom_res))


def populate_contacts_ofInterest(contacts, geneIntContacts, indexes):
    """Main function for creating the data frame of contacts for the genes of interest.

    contacts: A dictionary of all chromosome contacts.
    geneIntContacts: A data frame of empty gene contacts.
    indexes: A list og gene position tuples.
    """
    # geneIntContacts = geneIntContacts[~geneIntContacts.index.duplicated()]
    # geneIntContacts = geneIntContacts.loc[:,~geneIntContacts.columns.duplicated()]
    for i in range(len(indexes)):
        chrom1 = geneIntContacts.loc[indexes[i], "chr"].values[0]  # WTF is pandas .loc returning!!!!!!!
        # d = (chrom1, indexes[i][0], indexes[i][1])
        for j in range(i, len(indexes)):
            chrom2 = geneIntContacts.loc[indexes[j], "chr"].values[0]
            # r = (chrom2, indexes[j][0], indexes[j][1])
            dfCx = contacts[chrom1][chrom2]
            # Main condition that checkes in the flattened Hi-C map.
            if (indexes[i][1], indexes[j][1]) in dfCx.index:
                vc = dfCx.loc[(indexes[i][1], indexes[j][1])].values[0]
                if math.isnan(vc):
                    vc = 0.0
                # Replace the value in the large data frame.
                geneIntContacts.at[(indexes[i][0], indexes[i][1]), (indexes[j][0], indexes[j][1])] = vc



# TODO fix the argument ranges of accepted values from straw.
# TODO arguments: Add argument for organism.
parser = argparse.ArgumentParser(
    prog="pyna_collada",
    description="Tool to visualise Hi-C contacts from gene lists of interest.",
    epilog="Authors: Costas Bouyioukos, 2019-2022, Universite Paris Cite and UMR7216.")
parser.add_argument(
    "hic",
    type=str,
    metavar="hic_file",
    help="Filename (or path) of .hic file (NO option for STDIN)")
parser.add_argument(
    "genesAnnot",
    type=argparse.FileType("r"),
    default=None,
    metavar="gene_annot_coords_file",
    help="Filename (or path) of annotation file (usually from ENSEMBL) conaining at least the gene name, gene start and end, chromosome and strand")
parser.add_argument(
    "outfile",
    type=str,
    metavar="out_file",
    help="Filename (or path) of the results .html figure file (NO option for STDOUT)")
parser.add_argument(
    "-b",
    "--bin-size",
    help="Seelction of the bin size of the hi-c map (i.e. resolution). Default=10000",
    type=int,
    default=10000,
    dest="binSize",
    metavar="bin_size")
parser.add_argument(
    "-c",
    "--chromosomes",
    nargs="+",
    default="ALL",
    help="The chromosome names(s) of which we want to extract contacts. Deafult: ALL",
    metavar="chr_number",
    dest="chr")
parser.add_argument(
    "-d",
    "--downstream-offset",
    type=int,
    default=500,
    help="The number of bps to offset down-stream of TSS for defining the overlaping regions. Default: 500",
    metavar="downstream_off",
    dest="offsetD")
parser.add_argument(
    "-u",
    "--upstream-offset",
    type=int,
    default=1000,
    help="The number of bps to offset up-stream of TSS for defining the overlaping regions. Default: 1000",
    metavar="upstream_off",
    dest="offsetU")
parser.add_argument(
    "-i",
    "--interChromosomal",
    action="store_true",
    dest="inter",
    help="Flag to turn on inter- and intra-chromosomal contacts. Default: intra-chromosomal only")
parser.add_argument(
    "-n",
    "--normalisation",
    nargs="?",
    default="NONE",
    metavar="norm_meth",
    type=str,
    help="Choise of a normalisation method from the Juice-tools or hic-straw (One of VC, VC_SQRT, KR, Default: NONE)",
    dest="norm")
parser.add_argument(
    "-t",
    "--type",
    nargs="?",
    default="BP",
    metavar="Type",
    type=str,
    help="Choise of the hic-straw extracted measure. Default: BP",
    dest="type")
parser.add_argument(
    "-v",
    "--version",
    action="version",
    version="%(prog)s  v. {version}".format(version=__version__))


# Parse the command line arguments.
optArgs = parser.parse_args()
exit()
if optArgs.chr == "ALL":
    chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
else:
    chromosomes = optArgs.chr

# Compute the genes/contacts data frame
contacts = {}
for chrA in chromosomes:
    contacts[chrA] = extract_contacts(optArgs, chrA, chromosomes)
# Read and parse the features coordinate file.
with optArgs.genesAnnot as fh:
    next(fh)  # If the file contains a header!
    for l in fh:
        fields = l.split()
        if (fields[2] not in chromosomes):  #!!! Here we assume that the chromosome name is on the third column.
            continue
        # CAREFULL re-orienting genes to facilitate the analysis!!! we do not care so much *for the moment* for gene orientation.
        # if fields[4] < fields[3]:
        #    # Switch the start and end of a gene.
        #    tmp = fields[3]
        #    fields[3] = fields[4]
        #    fields[4] = tmp
        coords.append((fields[1], fields[2], int(fields[3]), int(fields[4])))  #!!! Here we assume the following column names ENSEMBL, GeneName, Chrom, GeneStart, GeneEnd OBLIGATORY!
labels = ["name", "chr", "start", "stop"]
# The genes of interest coordinates data frame
geneCoords = pd.DataFrame.from_records(gCoords, columns=labels)
# Sort the data frame according to chromosome and gene start site
geneCoords.sort_values(["chr", "start"], ascending=[True, True], inplace=True)
geneCoords.reset_index(drop=True, inplace=True)
# Find bins that overlap genes
intervals = []
for i, row in geneCoords.iterrows():
    startInt = int(row["start"]) // optArgs.binSize
    stopInt = int(row["stop"]) // optArgs.binSize
    interval = tuple([x * optArgs.binSize for x in range(startInt, stopInt + 1)])
    intervals.append(interval)
# Append the intervals into the data frame
geneCoords["intervals"] = intervals
# Expand the intervals / coordinates data frame
multiIntervs = []
names = []
bins = []
for i, row in geneCoords.iterrows():
    for j in row["intervals"]:
        if (row["name"] not in names) and (j not in bins):
            names.append(row["name"])
            bins.append(j)
            multiIntervs.append([row["name"], row["start"], row["stop"], row["chr"], j])
labels = ["name", "start", "stop", "chr", "bin"]
geneCoords = pd.DataFrame.from_records(multiIntervs, columns=labels)
# Prebuild the data frame of the matrix of genes of interest
# zip the name-X-bin columns to create the index tuples for rows and columns
indexes = list(zip(geneCoords["name"], geneCoords["bin"]))
# Build an empty data frame
geneIntContacts = pd.DataFrame(0, index=pd.MultiIndex.from_tuples(indexes), columns=pd.MultiIndex.from_tuples(indexes))
geneIntContacts.insert(0, "stop", list(geneCoords["stop"]))
geneIntContacts.insert(0, "start", list(geneCoords["start"]))
geneIntContacts.insert(0, "chr", list(geneCoords["chr"]))
# Main function to populate the data frame!
populate_contacts_ofInterest(contacts, geneIntContacts, indexes)
# FIXME Check if we really need logs or not, for the moment we use!
mm = np.log2(geneIntContacts.iloc[:, 3:].replace(0, np.nan))
mm = geneIntContacts.iloc[:, 3:].replace(0, np.nan)
# mm = geneIntContacts.iloc[:,3:]
# mm = mm.replace(0, np.nan)
indexes2 = ["{}-{}".format(a, b) for a, b in zip(geneCoords["name"], geneCoords["bin"])]
mm.index = indexes2
mm.columns = indexes2

# Ploting with plotly
# Transform pandas data frame to dictionary for the plotly visualisation.
ddMM = go.Heatmap(z=mm.to_numpy().tolist(), x=mm.index, y=mm.columns, hoverongaps=False)
fig = go.Figure(data=ddMM)
fig.update_layout(width=1600, height=1600,
                  title=f"Selection's contacts from chromosome {str(optArgs.chr[0])}",
                  xaxis_title="Gene names/Coordinates", yaxis_title="Gene names coordinates",
                  #legend_title="Normalised Contacts",
                  font=dict(family="Courier New, monospace",
                            size=14) # color="RebeccaPurple")
                  )
fig.update_yaxes(automargin=True)
fig.update_xaxes(automargin=True)
plotly.offline.plot(fig, filename=optArgs.outfile, auto_open=False)



# Ploting SNS DEPRICATED
# import matplotlib.pyplot as plt
# import seaborn as sns
# Setup
# sns.set(style="white")
# f, ax = plt.subplots(figsize=(13, 11))
# #cmap = sns.diverging_palette(220, 10, as_cmap=True)
# sns.color_palette(palette="OrRd")
# sns.heatmap(mm)
# mx = max(max(ax.get_ylim()), max(ax.get_xlim()))
# ax.set_ylim(mx, 0)
# ax.set_xlim(0, mx)
# plt.savefig(optArgs.outfile)
