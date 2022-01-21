from turtle import st
from dna_features_viewer import GraphicFeature, GraphicRecord
import pysam

chromosome = 'AGTU001.0001.c02'
start = 100000
end = 300000
bam = 'test_data/aligned.sorted.bam'
a = pysam.AlignmentFile(bam)
c = pysamstats.load_coverage(bam, chrom='Pf3D7_01_v3', pad=True)
reads = []

for read in a.fetch(chromosome, start, end):
    reads.append(read)

features = []
for read in reads:
    if read.is_reverse:
        strand = 1
    else:
        strand = -1
    gf = GraphicFeature(start=read.reference_start, end=read.reference_end, strand=strand,
                        color="#ffd700", label=read.query_name)
    features.append(gf)
record = GraphicRecord(sequence_length=end, features=features)
record = record.crop((start, end))

record.plot_on_multiple_pages("chunks.pdf",
                              nucl_per_line=end-start,
                              lines_per_page=10,
                              plot_sequence=False
                              )
