from dna_features_viewer import GraphicFeature, GraphicRecord
import pysam


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

    record.plot_on_multiple_pages(out,
                                  nucl_per_line=end-start,
                                  lines_per_page=10,
                                  plot_sequence=False
                                  )
