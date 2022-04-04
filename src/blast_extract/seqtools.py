"""
Nucleotide sequence related utilities.
"""

import typing as t

# TODO: replace with minimum viable self-written code,
#  and make Biopython optional dependency - used when available.
# From Biopython >>>
'''
BSD 3-Clause License
Copyright (c) 1999-2021, The Biopython Contributors All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following 
   disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products 
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''


def _maketrans(complement_mapping):
    """Make a python string translation table (PRIVATE).

    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.

    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.

    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    keys = "".join(complement_mapping.keys()).encode("ASCII")
    values = "".join(complement_mapping.values()).encode("ASCII")
    return bytes.maketrans(keys + keys.lower(), values + values.lower())


_ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
}
_dna_complement_table = _maketrans(_ambiguous_dna_complement)


# <<< End of Biopython

class Sequence:

    def __init__(self, name: str, fwdseq: str):
        self.name = name
        self.fwdseq = fwdseq.upper()
        self.revseq = reverse_complement(self.fwdseq)
        self.length = len(self.fwdseq)


def make_contigs_fasta(contigs: t.Iterable['Sequence']) -> str:
    result = ''
    for contig in contigs:
        result += f">{contig.name}\n{contig.fwdseq}\n"
    return result


def parse_genome_contigs(genome: t.TextIO, fsep) -> dict[str, 'Sequence']:
    """Parse a genome fasta file into a mapping of contig names and sequences."""
    contig_mapping = dict()
    for header, seq in parse_fasta(genome):
        name = header.split(fsep, 1)[0]  # First value in header assumed to be contig name
        contig_mapping[name] = Sequence(name=name, fwdseq=seq)
    return contig_mapping


def parse_fasta(fh: t.TextIO) -> t.Generator[t.Tuple[str, str], None, None]:
    """Parse a fasta file to tuples of header and sequence."""
    header = None
    seq = ""
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            # New header: output previous sequence if it exists, then clear sequence and set new header
            if seq:
                yield header, seq.replace(' ', '')
            seq = ""
            header = line[1:]
        elif header is not None:
            seq += line
        else:
            raise AssertionError(f"File did not start with FASTA header: {line}")
    # Output the last sequence
    if seq:
        yield header, seq


def complement(seq: str) -> str:
    data = seq.encode('ASCII')
    data = data.translate(_dna_complement_table)
    return data.decode('ASCII')


def reverse_complement(seq: str) -> str:
    return complement(seq=seq)[::-1]


#
#
# END OF FILE
