# Conflict And Duplication Identifier (CAnDI)

Scripts to investigate gene tree conflict with gene families.

**Python 3 version:** The program has also been converted to python3 and can be found in the folder Py3version. Although this version has been tested, there still may be some issues that arose from converting the program, and the python2 version used for the analyses in the paper is still the main one in the scripts folder for the time being.

## Introduction

This repository contains the CAnDI.py program designed to detect conflict using orthologs (orthologous gene trees) or gene families (Homologous gene trees). The program can either map a series of gene trees onto the species tree to identify the amount of conflict the gene trees have, map a species onto a gene tree to investigate where the gene tree has conflict, or search for patterns of conflict. A more detailed description of all the program's features may be found in the [manual](https://github.com/HollyMaeRobertson/CAnDI/blob/master/Manual_6-7.pdf).

License: GPL https://www.gnu.org/licenses/gpl-3.0.html

## Quick Start

For a simple analysis, using the prefix test1 for all output files (outfile prefix defaults to “out”), use the following code:

```python2 scripts/CAnDI.py --mode n --species_tree path/to/species_tree.tre --gene_folder path/to/gene_folder/ --outfile_prefix test1```

This generates a set of outfiles with the prefix "test1" that can be used to investigate conflict between the provided gene and species trees. For a tree with the nodes labelled with concordances, look for "test1_concord.tre" and the conflicts, "test1_conflict.tre". 

To quickly generate pie charts from the same data, using the same prefix of test1:

```python3 Pie.py -f . -p test1 -t test_total_analyzed.tre -o le -a 3```

To investigate a single homolog tree for concordance or conflict with the species tree, use reverse mode:

```python2 scripts/CAnDI.py --mode r --species_tree path/to/species_tree.tre --gene_tree path/to/gene_tree.tre```

Reverse mode only generates one outfile, called out_concon.tre (or equivalent prefix in place of “out”). This is a gene tree in newick format where each node is labelled as either in conflict (X), concordant (*), a duplication (D), or uninformative (U) with respect to the species tree.

For more information, run
```python2 scripts/CAnDI.py -h``` 
or 
```python3 OutputSummarizer/Pie.py -h```

Or, read the manual! 

In the event you run into any bugs or are for any other reason having trouble using CAnDI after reading the manual, do not hesitate to reach out to me at email hmr42@cam.ac.uk.  

## Publications
CAnDI is now published! See the paper here: https://doi.org/10.1093/sysbio/syaf028. 

In addition, CAnDI has been used to conduct conflict analysis in [this paper](https://doi.org/10.1093/molbev/msac044) about transcription factor evolution by Wheeler et al., 2022.  
