# Welcome to StrucFollower

## Introduction

This package helps you to identify deleted or inserted sequences in evolved bacterial strains.  
For example if you evolved antibiotic resistance to a bacterial strain this package helps you to identify sequences which were integrated or deleted in the genome during evolution. With it's support for GenBank files StrucFollower also tells you which products were inserted or deleted. StrucFollower optionally generates visualizations of alignments of inserted or deleted sequences which can help in assessing the confidence of the detected or inserted sequence.  

Additionally if you evolved your bacterial strain in co-culture with other bacterial strains, StrucFollower has a feature which identifies sequences which were integrated in the genome from the other bacterial strains of the co-culture.

This package has very few dependencies and because of it's simple design it's fast and reliable. 

## Installation

This package requires `SAMtools>=1.11` and `minimap2` in your PATH. If haven't installed those dependencies already you can install them with `conda`:

```
conda install -c bioconda samtools
conda install -c bioconda minimap2 
```

StrucFollower itself can be installed with pip:
```
pip install git+https://github.com/nahanoo/deletion_detection.git
```

StrucFolower creates three console scripts which are independently callable. Click on the console script name for the detailed documentation of the according sub-module. 

* [`detect_deletions`](deletion_detection) - Detects deleted sequences.  
* [`detect_insertions`](insertion_detection) - Detects inserted sequences.  
* [`detect_hgts`](horizontal_gene_transfer) - Detects sequences integrated in the genome from other bacterial strains.

## Input data

StrucFollower was developed and tested with PacBio long-read sequencing data and it's recommended to use high-quality assemblies. Alternatively to PacBio assemblies, hybrid-assembled genomes with Nanopre and Illumina data should work as well.

## Principle

All three sub-modules have a similar underlying mechanism. As an input the GenBank or the FASTA file of the ancestral strain and the mutated strain are required. Depending on the sub-module, either the ancestral or the mutated genome is chunked into smaller sequences using a sliding window. Typically a window size of 500 base-pairs and start-shift of 100 base-pairs is used. Below you can see an example of chunked sequences aligned to the genome itself:

![chunked_sequences](chunks.png)

These chunks are then aligned either to the ancestral or the mutated genome. For a more detailed description of the alignment check the according sub-module description.
