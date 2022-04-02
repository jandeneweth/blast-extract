"""
Hash a FASTA-format file of alleles.
"""

import sys
import typing
import Bio.SeqIO
import hashlib
import argparse


# -- Public Methods --

def iter_alleles(fh: typing.TextIO, fsep: str, field_nr: int):
    for record in Bio.SeqIO.parse(fh, 'fasta'):
        header = record.description
        fields = header.split(fsep)
        locus = fields[field_nr - 1] if len(fields) >= field_nr else None
        sequence = str(record.seq).lower()
        yield locus, sequence


def get_hash_fn(algorithm: str, hash_trim: int | None):
    try:
        hasher = hashlib.new(algorithm)
    except ValueError:
        raise RuntimeError(
            f"Invalid or unavailable hash algorithm for this platform: {algorithm}. "
            f"\nAvailable algorithms: {', '.join(sorted(hashlib.algorithms_available))}"
        )

    def do_hash(data):
        h = hasher.copy()
        h.update(data)
        result = h.hexdigest()
        if hash_trim is not None:
            result = result[:hash_trim]
        return result

    return do_hash


def run(fh: typing.TextIO, out: typing.TextIO, algorithm: str = 'md5', hash_trim: int | None = None, fsep: str = ' ', field: int = 1):
    hash_fn = get_hash_fn(algorithm=algorithm, hash_trim=hash_trim)
    for locus, sequence in iter_alleles(fh=fh, fsep=fsep, field_nr=field):
        if locus is None:
            continue
        hexhash = hash_fn(sequence)
        out.write(f"{locus}\t{hexhash}\n")


# -- Run as Script --

def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--algorithm', '-a', default='md5', type=str, help="The hashing algorithm used")
    parser.add_argument('--hash_trim', '-t', default=None, type=int, help="If specified, the hexadecimal representation of the output hashes are trimmed to the first X characters")
    parser.add_argument('--fsep', '-s', default=' ', type=str, help="The character(s) separating fields in the fasta headers")
    parser.add_argument('--field', '-f', default=1, type=int, help="The field in which the locus name can be found, starting at 1")
    parser.add_argument('--out', '-o', type=argparse.FileType('w'), default=sys.stdout, help="The output destination, defaults to STDOUT")
    parser.add_argument('ALLELES', help='The FASTA format input of alleles, defaults to STDIN', type=argparse.FileType('r'), default=sys.stdin)
    return parser


def main(args: list[str] | None = None):
    parser = get_argparser()
    ns = parser.parse_args(args=args)
    run(
        fh=ns.ALLELES,
        out=ns.out,
        algorithm=ns.algorithm,
        hash_trim=ns.hash_trim,
        fsep=ns.fsep,
        field=ns.field,
    )


if __name__ == '__main__':
    main()


#
#
# END OF FILE
