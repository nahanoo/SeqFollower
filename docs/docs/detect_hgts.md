# Detect horizontal gene transfers

## Alignment and detection

To detect horizontal gene transfers (HGTs) StrucFollower chunks the genome of the mutant into smaller chunks. 
Those chunks are aligned to the references of the strains of the co-culture.
Every alignment is scanned. If a sequence is mapped it's either conserved or originates from the strain of the alignment. Because we store the contig name and the step ID in the sequence name we can trace back the position of the sequence in the genome of the mutant. To check if the sequence was inserted in the mutant the sequence is aligned to the ancestor of the mutant. If the sequence aligns it was already present at the start of the experiment. The sequence suspected to be a product of a HGT is annotated. If the sequence was not inserted but is conserved among the strains this is often visible in the annotation. Such sequences are often annotated as conserved ribosomal sequences. On the other hand inserted sequences as a product of HGTs are often annotated as transposons.

## Usage

```
usage: detect_hgts [-h] [--plot] mutant ancestor references out_dir

Detect HGTs in bacteria evolved in co-cultures.

positional arguments:
  mutant      genbank file of the mutant
  ancestor    Fasta file of the ancestor. 
              Prefix of fasta file will be interpreted as strain name.
  references  Folder containing all genomes in fasta format of co-cultured bacteria.
              Prefix of fasta file will be interpreted as strain name.
  out_dir     output directory

optional arguments:
  -h, --help  show this help message and exit
  --plot      plot alignments and annotations
```

## Outputs

* `hgts.tsv` - Outputs all positions with sequences mapping to the provided references.  
The column `origin` stores to which reference the sequence aligns.
* `hgts.annotated.tsv` - Additionally stores the annotation for sequences mapping to the references.
* `plots` - Optional
    * `alignments` - Visualizes alignments used to detect hgts.
    * `annotations` - Visualizes the annotations of the detected sequences.