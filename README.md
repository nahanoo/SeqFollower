# SeqFollower

Detect inserted or deleted sequences from assemblies of evolved bacterial strains.
If evolved in co-culture you can use SeqFollower also to identify sequences inserted by horizontal gene transfers.
Because of it's simple design it's very fast taking less than 2 minutes to run locally.

## Documentation

[https://nahanoo.github.io/SeqFollower/](https://nahanoo.github.io/SeqFollower/)

## Installation 

This package requires `SAMtools>=1.11` and `minimap2` in your PATH. If you haven't installed those dependencies already you can install them with `conda`:

```
conda install -c bioconda samtools
conda install -c bioconda minimap2 
```

SeqFollower itself can be installed with pip:
```
pip install git+https://github.com/nahanoo/SeqFollower.git
```

SeqFollower creates three console scripts which are independently callable:

* [`detect_deletions`](detect_deletions) - Detects deleted sequences.  
* [`detect_insertions`](detect_insertions) - Detects inserted sequences.  
* [`detect_hgts`](detect_hgts) - Detects horizontal gene transfers.