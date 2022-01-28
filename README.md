# insertion_detection

Detect large insertions using a sliding window mapping algorithm.

## Installation 

```
git clone https://github.com/nahanoo/insertion_detection.git
cd insertion_detection
pip install .
```

This tool requires `samtools` and `minimap2` to be in your path.

## Usage

`detect_insertions -h`
```
usage: detect_insertions [-h] [--plot] anceteral mutant out_dir

Detect insertions in evolved bacterial strains. Fast and simpple.

positional arguments:
  anceteral   fasta file of the ancesteral strain
  mutant      genbank file of the mutated strain.
  out_dir     output directory

optional arguments:
  -h, --help  show this help message and exit
  --plot      if this flag is added the alignment of every insertion is
              plotted.
```

## Output

The output are two tsv files with the chromosome, position and length of the insertion. `insertions.annotated.tsv` also has the product names found in the insertion sequence. The alignments are plotted as well as the annotation of the insertion.

## Background

The reference sequence is chunked into 50000 base-pair sequences using a sliding window algorithm. Those chunks are aligned to the query mutant sequence. Regions with no coverage are not present in the reference and must therefore be insertions. The positions of these regions are stored and annotated.usage: detect_insertions [-h] [--plot] anceteral mutant out_dir

Detect insertions in evolved bacterial strains. Fast and simpple.

positional arguments:
  anceteral   fasta file of the ancesteral strain
  mutant      genbank file of the mutated strain.
  out_dir     output directory

optional arguments:
  -h, --help  show this help message and exit
  --plot      if this flag is added the alignment of every insertion is
              plotted.usage: detect_insertions [-h] [--plot] anceteral mutant out_dir

Detect insertions in evolved bacterial strains. Fast and simpple.

positional arguments:
  anceteral   fasta file of the ancesteral strain
  mutant      genbank file of the mutated strain.
  out_dir     output directory

optional arguments:
  -h, --help  show this help message and exit
  --plot      if this flag is added the alignment of every insertion is
              plotted.
