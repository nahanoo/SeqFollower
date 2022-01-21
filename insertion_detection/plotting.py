from dna_features_viewer import GraphicFeature, GraphicRecord
import pysam
from os.path import join
import pandas as pd


def plot_alignment(bam, chromosome, position, out):
    a = pysam.AlignmentFile(bam)
    padding = 50000
    if position - padding < 0:
        start = 0
    else:
        start = position - padding

    ref_length = a.get_reference_length(chromosome)
    if position + padding > ref_length:
        end = ref_length
    else:
        end = position + padding

    reads = []
    for read in a.fetch(chromosome, start, end):
        if (not read.is_unmapped) & (not read.is_secondary) & (read.mapq == 60):
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
    target = join(out, 'plots', '.'.join(
                [chromosome, str(position), 'pdf']))
    record.plot_on_multiple_pages(target,
                                  nucl_per_line=end-start,
                                  lines_per_page=10,
                                  plot_sequence=False
                                  )

def plot_products(bam,out):
    a = pysam.AlignmentFile(bam)
    f = join(out,'insertions.annotated.tsv')
    df = pd.read_csv(f,sep='\t')
    locations = set([(c,p) for c,p in zip(df['chromosome'],df['position'])])
    for location in locations:
        
    for i,row in df.iterrows():
        c = row['chromosome']
        p = row['position']
        l = row['length']
        start = row['start']
        end = row['end']
        product = row['product']
        padding = 10000
        if p - padding < 0:
            start = 0
        else:
            start = p - padding

        ref_length = a.get_reference_length(c)
        if p + padding > ref_length:
            end = ref_length
        else:
            end = p + padding
        
        features = []
        gf = GraphicFeature(start=read.reference_start, end=read.reference_end, strand=strand,
                            color="#ffd700", label=read.query_name)
        features.append(gf)

