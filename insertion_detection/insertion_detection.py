from os.path import join, exists
from os import mkdir, remove
from io import StringIO
from subprocess import call, run as r, DEVNULL, STDOUT
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from .plotting import plot_alignment, plot_genbank


class Insertion():
    def __init__(self, args):
        # Parsed files
        self.out_dir = args.out_dir
        self.reference = args.ancestral
        self.mutant_gbk = args.mutant

        # Step and window size for chunking
        self.step = 100
        self.window = 500

        # Mutant contigs and features
        self.mutant_contigs = [contig for contig in SeqIO.parse(
            self.mutant_gbk, 'genbank')]
        self.mutant_features = self.parse_genbank()
        self.mutant_fasta = join(self.out_dir, 'mutant.fasta')
        with open(join(self.out_dir, 'mutant.fasta'), 'w') as handle:
            SeqIO.write(self.mutant_contigs, handle, 'fasta')

        # Reference contigs
        self.reference_contigs = [contig for contig in SeqIO.parse(
            self.reference, 'fasta')]
        
        # Chunked mutant
        self.chunks = join(self.out_dir,
                      "chunked_sequences.fasta")

        # Mutant reference alignments
        self.bam = join(self.out_dir, "aligned.sorted.bam")

        # Out dataframes
        self. insertions = pd.DataFrame(columns=['chromosome', 'position', 'length'])
        self.annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'product'])

        # Items to trash
        self.trash = []
        self.trash.append(self.mutant_fasta)

        if not exists(join(self.out_dir, 'plots')):
            mkdir(join(self.out_dir, 'plots'))

    def parse_genbank(self):
        genbank = {contig.id: {} for contig in self.mutant_contigs}
        for contig in self.mutant_contigs:
            for feature in contig.features:
                try:
                    start = feature.location.start
                    end = feature.location.end
                    product = feature.qualifiers['product']
                    genbank[contig.id][(start, end)] = product[0]
                except KeyError:
                    pass
        return genbank

    def chunker(self, seq, window_size, step):
        """Creates chunks of a sequence. window_size defines
        chunk size and step the amount of basepairs the window
        is moved forward."""
        # List which stores all chunks
        seqs = []
        seqlen = len(seq)
        self.step = step
        for counter, q in enumerate(range(0, seqlen, step)):
            # Returns ether entire sequence or window depending on sequence length
            j = seqlen if q + window_size > seqlen else q + window_size
            chunk = seq[q:j]
            # Add chunk id to sequence id
            chunk.id = chunk.id + "." + str(counter)
            seqs.append(chunk)
            if j == seqlen:
                break
        return seqs

    def chunk_assembly(self):
        """Chunks an assembly of multiple contigs into different 
        chunks using a sliding window algorightm (see chunker function)."""
        assembly_chunks = []
        for contig in self.reference_contigs:
            # Creates chunks of every contig
            assembly_chunks += self.chunker(contig, self.window, self.step)
        # Dumps chunks to fasta
        with open(self.chunks, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")
        self.trash.append(self.chunks)

    def mapper(self):
        """Maps chunked sequence to all ancesteral genomes
        experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        reads = join(self.out_dir, "chunked_sequences.fasta")
        sam = join(self.out_dir, "aligned.sam")
        cmd = [
            "minimap2",
            "-ax",
            "asm5",
            self.mutant_fasta,
            reads,
            ">",
            sam,
        ]
        # Calling minimap and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'sort', '-o', self.bam, sam]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'index', self.bam]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        if sam not in self.trash:
            self.trash.append(sam)
            self.trash.append(self.bam)
            self.trash.append(self.bam+'.bai')
        

    def get_insertions(self):
        """Get all areas with no coverage."""
        cmd = ['samtools', 'depth', '-aa', '-J', '-Q', '0', self.bam]
        process = r(cmd, capture_output=True)
        df = pd.read_csv(StringIO(process.stdout.decode()), sep='\t')
        df.columns = ['chromosome', 'position', 'coverage']
        df = df[df['coverage'] == 0]
        self.insertions = self.concat_insertions(df)
        # Switching to zero based output
        self.insertions['position'] = self.insertions['position'] - 1

    def concat_insertions(self, df):
        """Concatting seperated positions to insertions with n length."""
        i = -1
        prev_pos = 0
        prev_contig = None
        for contig, pos in zip(df['chromosome'], df['position']):
            if (prev_contig == contig) & (pos - 1 == prev_pos):
                self.insertions.at[i, 'length'] += 1
            else:
                i += 1
                self.insertions.loc[i] = [contig, pos, 1]
            prev_pos = pos
            prev_contig = contig
        return self.insertions

    def annotate(self):
        """Annotated insertions."""
        i = 0
        for counter, row in self.insertions.iterrows():
            c = row['chromosome']
            p = row['position']
            l = row['length']
            products = self.annotate_position(c, p, l)
            for product in products:
                self.annotated.loc[i] = [c, p, l, product]
                i += 1
        self.annotated = self.annotated.drop_duplicates()

    def annotate_position(self, c, p, l):
        """Returns products in a region in a genbank."""
        products = []
        for (start, end), product in self.mutant_features[c].items():
            if not set(range(start, end)).isdisjoint(range(p, p+l)):
                products.append(product)
        return products

    def dump(self):
        """Write insertions to tsv."""
        self.insertions.to_csv(
            join(self.out_dir, 'insertions.tsv'), sep='\t', index=False)
        self.annotated.to_csv(
            join(self.out_dir, 'insertions.annotated.tsv'), sep='\t', index=False)

    def plot_insertions(self):
        """Plot alignments."""
        out = join(self.out_dir, 'plots', 'alignments')
        if not exists(out):
            mkdir(out)

        for chromosome, position, length in zip(self.insertions['chromosome'], self.insertions['position'],
                                                self.insertions['length']):
            plot_alignment(self.bam, chromosome, position,
                           length, self.window, out)

    def plot_annotation(self):
        """Plot inserted products."""
        out = join(self.out_dir, 'plots', 'annotations')
        if not exists(out):
            mkdir(out)
        for i, row in self.insertions.iterrows():
            plot_genbank(
                self.mutant_contigs, row['chromosome'], row['position'],
                row['position']+row['length'], out)

    def clean(self):
        """Clean temporary files."""
        for item in self.trash:
            remove(item)
