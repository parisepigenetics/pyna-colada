"""Modle to facilitate hic contacts extraction and visualisation.
"""

import math
import numpy as np
import pandas as pd
import multiprocessing as mp
# TODO use functools for the expensive computation parts
#import functools
from hicstraw import straw
import plotly.graph_objects as go

import time    # for timing
import pprint  # for testing

def timing(f):
    """Wrapper to time functions.py

    Works as a decorator. Taken from https://stackoverflow.com/questions/5478351/python-time-measure-function
    """
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print("{:s} function took {:.3f} ms".format(f.__name__, (time2 - time1) * 1000.0))
        return ret
    return wrap


def parse_annotation(optArgs):
    """Extract annotation, compute gene start site and overalp intervals

    :param annotFile: An ENSEMBL annotation file with (at least) gene name, gene start and end sites and strand
    :type annotFile: An open ready to read filehandler
    :return: A pandas data frame with the ENSEMBL ids as row names and gene common name, begin, end of interaction sites as columns
    :rtype: pandas.DataFrame
    """
    # TODO for future release... fetch these annotations from ENSEMBL online API
    # Read and parse the features coordinate file
    ga = pd.read_csv(optArgs.genesAnnot, delimiter='\t', low_memory=False)
    # Setting up the index
    #idx = list(ga.loc[:,"Gene_stable_ID"].values)
    #ga.index = idx  # No need to change the index perhaps TODO re-check it later in the production
    #ga = ga.drop("Gene_stable_ID", axis=1)
    # TODO do some post-checking for NaNs in strad, chromosome and gene start coords
    # Calculate the overlaping intervals
    # Get the positive strand indexes
    ipos = ga[ga["Strand"]==1].index.tolist()
    # Getting the negative strand indexes2
    ineg = ga[ga["Strand"]==-1].index.tolist()
    # Calculate the overlaping intervals based on gene start strand and the up- and down- offsets
    ga.loc[ipos, "InterStart"] = ga.loc[ipos, "Gene_start"] - optArgs.offsetU
    ga.loc[ipos, "InterEnd"] = ga.loc[ipos, "Gene_start"] + optArgs.offsetD
    ga.loc[ineg, "InterStart"] = ga.loc[ineg, "Gene_end"] - optArgs.offsetD
    ga.loc[ineg, "InterEnd"] = ga.loc[ineg, "Gene_end"] + optArgs.offsetU
    # Calculate the corresponding interacting bins
    ga.loc[:,"BinStart"] = round(ga.loc[:,"InterStart"] / optArgs.binSize) * optArgs.binSize
    ga.loc[:,"BinEnd"] = round(ga.loc[:,"InterEnd"] / optArgs.binSize) * optArgs.binSize
    # Covert the boundaries to integers
    ga = ga.astype({"InterStart" : int, "InterEnd" : int, "BinStart" : int, "BinEnd" : int})
    # Sort the data frame according to chromosome and interaction region start "InterStart"
    ga.sort_values(["Chromosome", "InterStart"], ascending=[True, True], inplace=True)
    ga.reset_index(inplace=True, drop=True)
    return ga


def expand_gene_annotation_bins(gadf, optArgs):
    """Expandds the multiple intervals per gene or the multiple genes per interval.

    Takes the genee annotation data frame from parse_annotation and expands the intervals such that the many-to-many relations between genes and HiC bins become one-to-one (with mupliple genes per bin if it is necessary).
    :param gadf: A gene annotation pandas data frame which contains the interaction region of interest for each gene as columns names "InterStart" and "InterEnd".
    :type gadf: pandas.DataFrame
    :return: A pandas Dataframe with the intervals (i.e. HiC file bins) for each gene, multi-indexed by chromosome-/-gene.
    :rtype: pandas.DataFrame
    """
    # Get the chromosomes
    chroms_DFs = []
    chroms = set(gadf["Chromosome"].values.tolist())
    for chr in chroms:
        # Generate an empty sub-DataFrame for the chromosome
        chrDF = pd.DataFrame()
        # Get the chomosome slice form the original DataFrame
        gaSlice = gadf.loc[gadf["Chromosome"]==chr, :]
        # Iterrate over rows
        for i, row in gaSlice.iterrows():
            # Set the condition to generate the expanded rows per gene over an interval
            binTemp = row["BinStart"]
            while binTemp <= row["BinEnd"]:  # This while loop produces the expansion!
                adds = pd.Series({"BinPlot" : binTemp, "GenePlot" : row["Gene_name"]})
                newRow = pd.concat([row, adds])
                chrDF = pd.concat([chrDF, newRow.to_frame().T])
                binTemp = binTemp + optArgs.binSize
        chroms_DFs.append(chrDF)
    # Collect all the chromosomes DFs in one
    gaexp = pd.concat(chroms_DFs)
    # Sort it out and reset the index
    gaexp.sort_values(["Chromosome", "InterStart"], ascending=[True, True], inplace=True)
    gaexp.reset_index(inplace=True, drop=True)
    return gaexp


@timing
def get_contacts_frame(optArgs, chrA, chrB):
    """Extract the contact matrix from a .hic fle using the straw interface.

    :return: A double indexed (chromosomes-/-bins) data frame with the contacts.
    :rtype: pandas.DataFrame
    """
    res = straw("observed", optArgs.norm, optArgs.hicFile, chrA, chrB, optArgs.type, optArgs.binSize)
    # The new API of hicstraw returns a list of contactRecord objects so we extract the bins and counts by list comprehension
    data = [(r.binX, r.binY, r.counts) for r in res]
    cont = "counts_" + chrA + "x" + chrB + "_" + "Norm" + optArgs.norm
    dc = pd.DataFrame(data=data, columns=["binX", "binY", cont])
    dc = dc.set_index(["binX", "binY"])
    return dc


def extract_contacts(optArgs, chrA, chromosomes):
    """The paralelisation wrapper of the contacts extract process.Η μπρο­σού­ρα της Αλε­ξάν­δρα Κολ­λο­ντάι (1872–1952) «Κομ­μου­νι­σμός και Οικο­γέ­νεια» πρω­το­δη­μο­σιεύ­θη­κε το 1920. Το κεί­με­νο που ακο­λου­θεί, μαζί με την Εισα­γω­γή, είναι από την ομώ­νυ­μη μπρο­σού­ρα των εκδό­σε­ων «Λάβα» του 1974, που έχει εξαντληθεί.



    Works one chromosome at a time: i.e. it paralelises all the contacts beteen all the chromosome pairs in the user defined list of input.
    """
    # NOT IN USE AS THE TOOL WORKS WITH A SINGLE CHROMOSOME AT A TIME FOR THE moment
    # TODO to be developed as a prrallelisation function IFF the tool is expanded to visualise many cromosomes.
    pass

@timing
def populate_contacts_ofInterest(contacts, gaExpDf):
    """Main function for creating the data frame of contacts for the genes of interest.

    :param contacts: A dictionary of chromosome of chromosomes contacts.
    :type contacts: dict
    :param gaExpDf: The expanded gene annotation data frame.
    :type gaExpDf: pandas.DataFrame
    :return: A list of pandas DataFrames with the contacts of interest.
    :rtype: list of pandas.DataFrame
    """
    contacts_OfI_list = []
    # Get the annotation ines corresponding to the chromosomes
    for chr1 in contacts.keys():
        for chr2 in contacts[chr1].keys():
            if chr1 == chr2:
                lines = gaExpDf.loc[gaExpDf["Chromosome"]==chr1, :]
            else:  # The case of inter-chromsomal contacts
                lines = gaExpDf.loc[(gaExpDf["Chromosome"]==chr1) | (gaExpDf["Chromosome"]==chr2), :]
            lines.reset_index(inplace=True, drop=True)
            # Get the corresponding contacts data
            dfCx = contacts[chr1][chr2]
            # Generate an emty data frme
            geneIntContacts = pd.DataFrame(0, index=lines["BinPlot"].values.tolist(), columns=lines["BinPlot"].values.tolist())  # FIXME this need to have a double index for inter-chromosomal contacts.
            # Populate the geneIntContacts DataFrame
            for i in range(len(lines)):
                for j in range(i, len(lines)):
                    binx = lines["BinPlot"][i]
                    biny = lines["BinPlot"][j]
                    # Find the counts for the bins of interest
                    if (binx, biny) in dfCx.index:
                        vc = dfCx.loc[(binx, biny)].values[0]
                        if math.isnan(vc):
                            vc = 0.0
                        # Replace the value in the large data frame.
                        geneIntContacts.at[binx, biny] = vc
            geneIntContacts = geneIntContacts.replace(0, np.nan)
            newIndexes = [f"{a}-{int(b/1000):,}Kb" for a, b in zip(lines["Gene_name"], lines["BinPlot"])]
            geneIntContacts.index = newIndexes
            geneIntContacts.columns = newIndexes
            contacts_OfI_list.append(geneIntContacts)
    return contacts_OfI_list


def matrx_trace(moiList, l=0):
    """Produce the plottly trace from a list of the matrix of contacts of interest.

    """
    trace = go.Heatmap(z=np.array(moiList[l]), x=moiList[l].columns, y=moiList[l].index,
                    hoverongaps=False, colorscale="bluered")
    return trace


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
