"""
Nucleotide sequence related utilities.
"""

import typing as t

import Bio.Seq
import Bio.SeqIO


class Contig:

    def __init__(self, name: str, fwdseq: str):
        self.name = name
        self.fwdseq = fwdseq.upper()
        self.revseq = reverse_complement(self.fwdseq)
        self.length = len(self.fwdseq)


def make_contigs_fasta(contigs: t.Iterable['Contig']) -> str:
    result = ''
    for contig in contigs:
        result += f">{contig.name}\n{contig.fwdseq}\n"
    return result


def parse_genome_contigs(genome: t.TextIO, fsep) -> dict[str, 'Contig']:
    """Parse a genome fasta file into a mapping of contig names and sequences."""
    contig_mapping = dict()
    for seqrecord in Bio.SeqIO.parse(genome, 'fasta'):
        name = seqrecord.description.split(fsep, 1)[0]  # First value in header assumed to be contig name
        seq = str(seqrecord.seq).lower()
        contig_mapping[name] = Contig(name=name, fwdseq=seq)
    return contig_mapping


def reverse_complement(seq: str) -> str:
    # TODO: reimplement this to remove biopython dependency
    return Bio.Seq.reverse_complement(seq)


#
#
# END OF FILE
