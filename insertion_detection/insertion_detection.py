from Bio import SeqIO
from os.path import join
from subprocess import call, run as r, DEVNULL, STDOUT
from io import StringIO
import pandas as pd
from plotting import plot_alignment


class Insertion():
    def __init__(self):
        self.out_dir = 'test_data'
        self.step = 10000
        self.window = 50000
        self.mutant_gbk = 'test_data/assembly.gbk'
        self.mutant_contigs = [contig for contig in SeqIO.parse(
            self.mutant_gbk, 'genbank')]
        self.mutant_fasta = join(self.out_dir, 'mutant.fasta')
        with open(join(self.out_dir, 'mutant.fasta'), 'w') as handle:
            SeqIO.write(self.mutant_contigs, handle, 'fasta')
        self.reference = 'test_data/reference.fasta'
        self.ref_contigs = [contig for contig in SeqIO.parse(
            self.reference, 'fasta')]
        self.insertions = None
        self.genbank = None
        self.bam = join(self.out_dir, "aligned.sorted.bam")
        self.annotated = pd.DataFrame(
            columns=['chromosome', 'position', 'length', 'start', 'end', 'protein'])

    def parse_genbank(self):
        self.genbank = {contig.id: {} for contig in self.mutant_contigs}
        for contig in self.mutant_contigs:
            for feature in contig.features:
                try:
                    start = feature.location.start
                    end = feature.location.end
                    product = feature.qualifiers['product']
                    self.genbank[contig.id][(start, end)] = product[0]
                except KeyError:
                    pass

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
        for contig in self.ref_contigs:
            # Creates chunks of every contig
            assembly_chunks += self.chunker(contig, self.window, self.step)
        target = join(self.out_dir,
                      "chunked_sequences.fasta")
        # Dumps chunks to fasta
        with open(target, "w") as handle:
            SeqIO.write(assembly_chunks, handle, "fasta")

    def mapper(self):
        """Maps chunked sequence to all ancesteral genomes
        experiment with minimap2.
        Minimap2 settings are set to accureate PacBio reads."""
        reads = join(self.out_dir, "chunked_sequences.fasta")
        sam = join(self.out_dir, "aligned.sam")
        cmd = [
            "minimap2",
            "--secondary=no",
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
        # Calling minimap and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)

    def get_insertions(self):
        cmd = ['samtools', 'depth', '-aa', '-J', '-Q', '60', self.bam]
        process = r(cmd, capture_output=True)
        df = pd.read_csv(StringIO(process.stdout.decode()), sep='\t')
        df.columns = ['chromosome', 'position', 'coverage']
        df = df[df['coverage'] == 0]
        self.insertions = self.concat_insertions(df)
        self.insertions['position'] = self.insertions['position'] - 1
        self.insertions.to_csv(join(self.out_dir, 'insertions.tsv'), sep='\t')

    def concat_insertions(self, df):
        out = pd.DataFrame(columns=['chromosome', 'position', 'length'])
        i = -1
        prev_pos = 0
        prev_contig = None
        for contig, pos in zip(df['chromosome'], df['position']):
            if (prev_contig == contig) & (pos - 1 == prev_pos):
                out.at[i, 'length'] += 1
            else:
                i += 1
                out.loc[i] = [contig, pos, 1]
            prev_pos = pos
            prev_contig = contig
        return out

    def annotate(self):
        self.parse_genbank()
        i = 0
        for counter, row in self.insertions.iterrows():
            c = row['chromosome']
            p = row['position']
            l = row['length']
            for (start, end), product in self.genbank[c].items():
                if (start >= p) & (end <= p+l):
                    self.annotated.loc[i] = [c, p, l, start, end, product]
                    i += 1
        self.annotated.to_csv(
            join(self.out_dir, 'insertions.annotaed.tsv'), sep='\t')

    def plot_insertions(self):
        for chromosome, position in zip(self.insertions['chromosome'], self.insertions['position']):
            target = join(self.out_dir,'.'.join([chromosome,str(position),'pdf']))
            plot_alignment(self.bam,chromosome,position,target)

i = Insertion()
# i.mapper()
no_cov = i.get_insertions()
i.annotate()
i.plot_insertions()
