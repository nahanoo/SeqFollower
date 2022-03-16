# Welcome to StrucFollower

## Introduction

This package helps you to identify deleted or inserted sequences in evolved bacterial strains.  
For example if you evolved antibiotic resistance to a bacterial strain this package helps you to identify sequences which were integrated or deleted in the genome during evolution. With its support for GenBank files StrucFollower also tells you which products were inserted or deleted. StrucFollower optionally generates visualizations of alignments of inserted or deleted sequences which can help in assessing the confidence of the detected sequence.  

Additionally if you evolved your bacterial strain with other bacterial strains in co-culture, StrucFollower has a feature which identifies sequences which were integrated in the genome from the other bacterial strains.

This package has very few dependencies and because of the simple design it's fast and reliable. 

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

* [`detect_deletions`](detect_deletions) - Detects deleted sequences.  
* [`detect_insertions`](detect_insertions) - Detects inserted sequences.  
* [`detect_hgts`](horizontal_gene_transfer) - Detects sequences integrated in the genome from other bacterial strains.

## Input data

StrucFollower was developed and tested with PacBio long-read sequencing data and it's recommended to use high-quality assemblies. As an alternative to PacBio assemblies, hybrid-assemblies with Nanopre and Illumina data should work as well.

## Principle

All three sub-modules have a similar underlying mechanism. As an input the GenBank or the FASTA file of the ancestral strain (ancestor) and the mutated strain (mutant) are required. Depending on the sub-module, either the genome of the ancestor or the mutant is chunked into smaller sequences using a sliding window. Typically a window size of 500 base-pairs and a start-shift of 100 base-pairs is used. Below you can see an example of chunked sequences aligned to the genome itself:

![chunked_sequences](chunks.png)

The sequence name consists of the contig name and an enumerated counter of the chunk. This helps us to verify deletions or insertions with the generated visualizations. Below you can see an example of a visualization for a detected insertion where the chunked genome of the mutant was aligned to the genome of the ancestor.

![insertion](insertion.png)

The gap in the alignment means that the sequence is only present in the mutant but not in the ancestor. Because of the continous counters visible in the sequence names we can see that the inserted sequence is flanked by high-quality alignments which increases the confidence in the detected insertion.
