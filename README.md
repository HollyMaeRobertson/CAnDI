# Conflict And Duplication Identifier (CAnDI)
Scripts to investigate gene tree conflict with gene families. 

**Update: CAnDI is now published!** Please cite https://doi.org/10.1093/sysbio/syaf028 where applicable :D. 

## Introduction
This repository contains the CAnDI.py program designed to detect conflict using orthologs (orthologous gene trees) or gene families (homologous gene trees). The program can map a series of gene trees onto the species tree to identify the amount of conflict the gene trees have, map a species onto a gene tree to investigate where the gene tree has conflict, or search for patterns of conflict. A more detailed description of all the program's features may be found in the [manual](https://github.com/HollyMaeRobertson/CAnDI/blob/master/Manual_6-7.pdf).

License: GPL https://www.gnu.org/licenses/gpl-3.0.html

## Quick Start
### Install and Dependencies
To install CAndI, just clone the repository into your chosen folder:

```git clone https://github.com/HollyMaeRobertson/CAnDI.git```

The scripts should run immediately, no further setup necessary! CAnDI can run on Python3 or Python2 (Python3 recommended). 

### Counting conflicts on the species tree
For a simple analysis, use the following code:

```python3 Py3version/CAnDI.py --mode n --species_tree path/to/species_tree.tre --gene_folder path/to/gene_folder/```

where *path/to/species_tree.tre* and *path/to/gene_folder/* are replaced with the paths to the species tree and the directory of gene trees being investigated respectively.

This generates a set of outfiles with the default prefix "out" that can be used to investigate conflict between the provided gene and species trees. For a newick tree with the nodes labelled with concordances look for *out_concord.tre*, and for the conflicts, *out_conflict.tre*. 

To quickly generate pie charts from the same data, using the custom prefix of test1:

```python3 Pie.py -f . -p test1 -t test_total_analyzed.tre -o le -a 3```.

These pie charts will be labelled with node numbers, which can be mapped onto the tree labelled with node numbers, *out_labels.tre*.

### Mapping conflicts to a single homolog tree
To investigate a single homolog tree for concordance or conflict with the species tree, use reverse mode:

```python3 Py3version/CAnDI.py --mode r --species_tree path/to/species_tree.tre --gene_tree path/to/gene_tree.tre```.

Reverse mode only generates one outfile, called *out_concon.tre* (or your chosen prefix in place of "*out*"). This is a gene tree in newick format where each node is labelled as either in conflict (X), concordant (*), a duplication (D), or uninformative (U) with respect to the species tree.

### Manual and support
For more information, run
```python3 Py3version/CAnDI.py -h``` 
or 
```python3 OutputSummarizer/Pie.py -h```.

Or, read the  [manual](https://github.com/HollyMaeRobertson/CAnDI/blob/master/Manual_6-7.pdf)! 

In the event you run into any bugs or are for any other reason having trouble using CAnDI after reading the manual, do not hesitate to reach out to me at hmr42@cam.ac.uk.  

## Publications
CAnDI is now published! See the paper here: https://doi.org/10.1093/sysbio/syaf028. 

In addition, CAnDI has been used to conduct conflict analysis in [this paper](https://doi.org/10.1093/molbev/msac044) about transcription factor evolution by Wheeler et al., 2022.  

## Python 3 version 
The program has been converted to python3 and can be found in the folder Py3version. Although this version has been tested, there still may be some issues that arose from converting the program, and the python2 version used for the analyses in the paper is still in the scripts folder. 
