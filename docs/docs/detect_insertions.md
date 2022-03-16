# Detect insertions

## Principle

To detect insertions StrucFollower chunks the genome of the ancestor into smaller chunks. These chunks are aligned to the genome of the mutant and the alignment is scanned for areas with no coverage. In areas with no coverage the given sequence is present in the genome of the mutant but absent the genome of the ancestor.

## Usage

```
usage: detect_insertions [-h] [--plot] ancestral mutant out_dir

Detect insertions in evolved bacterial strains.

positional arguments:
  ancestral   fasta file of the ancestral strain
  mutant      genbank file of the mutated strain
  out_dir     output directory

optional arguments:
  -h, --help  show this help message and exit
  --plot      plot alignments and annotations
```

Example:
```
detect_insertions --plot ancestral.fasta mutant.gbk ./
```

## Outputs

* `insertions.tsv` - Outputs all positions with inserted sequences
* `insertions.annotated` - Additionally outputs which products were inserted
* `plots` - Optional
    * `alignments` - Visualizes insertions by plotting the alignments of the sequence chunks of the ancestor aligned to the genome of the mutant. 
    * `annotations` - Visualizes the annotation of the inserted sequence.