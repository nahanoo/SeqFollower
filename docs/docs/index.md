# Welcome to StrucFollower

Identify structural changes in evolved bacterial strains.

## Introduction

This package helps you to identify deleted or inserted sequences in evolved bacterial strains.
For example if you have evolved antibiotic resistance to a bacterial strain this package helps you to identify sequences which were integrated or deleted during evolution. With it's support for GenBank files StrucFollower also tells you which products were inserted or deleted. To help you in assessing if the detected or inserted sequences are likely to be true StrucFollower outputs plots of alignments of regions with deleted or inserted sequences.  
Additionally if you evolved your bacteria in co-culture with other bacterial strains, StrucFollower also has a feature which identifies sequences which were integrated in the genome from the other bacterial strains in the co-culture.
This package has very few dependencies and because of it's simple design it's fast and reliable. 

## Installation

This package requires `SAMtools>=1.11` and `minimap2` in your PATH. If you don't have those dependencies installed already you can install them with `conda`:

```
conda install -c bioconda samtools
conda install -c bioconda minimap2 
```

StrucFollower can be installed with pip:
```
pip install git+https://github.com/nahanoo/deletion_detection.git
```

After installing with pip StrucFolower creates three console scripts which are independently callable:

* [`detect_deletions`](#detect_deletions) - Detects deleted sequences.  
* `detect_insertions` - Detect inserted sequences.  
* `hgt` - Detect sequences integrated in the genome from other bacterial strains.

## Principles

All three sub-modules have a similar underlying principle. As an input the GenBank or the FASTA file of the ancestral strain and the mutated strain are required. Depending on the sub-module, either the ancestral of the mutated genome is chunked into smaller sequences using a sliding window. Typically a window size of 500 base-pairs and start shift of 100 base-pairs is chosen. Below you can see an example chunked sequences aligned to the genome itself:
![chunked_sequences](chunks.png)

## detect_deletions

