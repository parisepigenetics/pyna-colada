# Pyna Colada
-------------

Pyna colada is a python module and application using 'straw' and 'juice-tools' to provide a more useful API to work with .hic Hi-C experiments files.
------------------
Developed by [Costas Bouyioukos](https://github.com/cbouyio) at Universite Paris Cite and UMR7216 [Paris Epigenetics](https://github.com/parisepigenetics) since Dec. 2019

## Description

Pyna colada includes a library (pyna_colada) providing easy access and visualisation to .hic files on top of the [hic-straw](https://github.com/aidenlab/straw) library and an executable script (pyna-collda) that visualises contacts between any gene (or regions) list of interest by simply using an ENSEMBL (or other) annotation file.

### Installation
Just run the `setup.py` script to install the module and the executable in your path of interest. Typically:

```python setup.py install --user```

### Usage
A typical run of pyna-collada can be performed by specifying two input files as positional arguments and an output .html file with the interactive plotly visualisation:

```pyna-collada <.hic_input_file> <annotation_file> <output_html_file> [options]```

Where:
  * <hic_input_file>: is a binary .hic file usually the outcome of HiC experiment.
  * <annotation_file>: is a table (tab delimited) file where each line represents a gen (or region) of interest and with 3 (at least) mandatory columns "Gene_ID" (any unique identifier), "Gene_start", "Gene_end" (the beginning and end of gene in genomic coordinates), "Chromosome" (the chromosome where the gene resides), "Strand" (the + or - strand), "Gene_name" (the common gene name for visualisation purposes)
  * <output_html_file>: the path/name of the output file.
 
A description of the available options (and their default values) is fully accessible by the `--help` or `-h` command line flags.

## Prerequisites

1. pandas
2. hic-straw
3. numpy
4. plotly

## Contact
Developer/maintainer Costas Bouyioukos https://github.com/cbouyio, any comments or question please drop a line at: costas.bouyioukos@u-paris.fr