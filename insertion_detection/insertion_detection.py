from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os.path import join
from os import listdir
from subprocess import call, DEVNULL, STDOUT


class Insertion():
    def __init__(self):
        self.step = 10000
        self.window = 50000
        mutant = 'test_data/assembly.gbk'
        self.reference = 'test_data/reference.fasta'
        self.contigs = [contig for contig in SeqIO.parse(
            mutant, 'genbank')]
        self.out_dir = 'test_data'

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
        for contig in self.contigs:
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
        bam = join(self.out_dir, "aligned.sorted.bam")
        cmd = [
            "minimap2",
            "-ax",
            "asm5",
            self.reference,
            reads,
            ">",
            sam,
        ]
        # Calling samtools and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools', 'sort', '-o', bam, sam]
        # Calling minimap and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)
        cmd = ['samtools','index',bam]
        # Calling minimap and surpressing stdout
        call(" ".join(cmd), shell=True, stdout=DEVNULL,
             stderr=STDOUT)


i = Insertion()
i.chunk_assembly()
i.mapper()
