# Insertion Detection

## Principle

To detect insertions StrucFollower chunks the genome of the ancestral strain into smaller chunks. These chunks are aligned to the genome of the mutated strain. The alignment is scanned for areas with no coverage. In areas with no coverage the given sequence is present in the genome of the mutated strain but absent the genome of the ancestral strain.

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
  --plot      if this flag is added the alignment of every insertion is plotted
```

Example:
```
detect_insertions --plot ancestral.fasta mutant.gbk ./
```

## Outputs

* `insertions.tsv` - Outputs all positions with inserted sequences
* `insertions.annotated` - Additionally outputs which products were inserted
* `plots` - Optional
    * `alignments` - Visualizes alignments of areas with no coverage. 
    * `annotations` - Visualizes the annotation of the inserted sequence.

Visualization of a detected insertion
![insertion](insertion.png)
The base-track is the genome of the mutated strain and the aligned sequences the chunks of the genome of the ancestral strain.
Because we aligned the genome of the ancestral strain to the genome of the mutated strain, insertions are visible as areas with no coverage.  
The ID of the chunk consists of the contig name and the enumerated counter of the chunk. As you you can see the highest counter on the right side of the gap is 14429. The counter is continued on the left side with 14429. This continuity ensures us that the gap in the alignments is due to a inserted sequence. In this example we can say pretty confidently that we identified an inserted sequence. If the counters wouldn't be continuous or even come from different contigs the confidence would be lower.