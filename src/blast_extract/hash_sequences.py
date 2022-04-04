"""
Hash a FASTA-format file of alleles.
"""

import sys
import hashlib
import argparse
import typing as t

from blast_extract import seqtools


# -- Public Methods --

def iter_sequences(fh: t.TextIO, fsep: str, field_nr: int):
    for header, sequence in seqtools.parse_fasta(fh):
        fields = header.split(fsep)
        locus = fields[field_nr - 1] if len(fields) >= field_nr else None
        sequence = str(sequence)
        yield locus, sequence


def get_hash_fn(algorithm: str, hash_trim: int | None) -> t.Callable[[bytes], str]:
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


def run(
        fh: t.TextIO,
        out: t.TextIO,
        algorithm: str = 'md5',
        hash_trim: int | None = None,
        algprefix: bool = False,
        validate: bool = True,
        fsep: str = ' ',
        field: int = 1
):
    hash_fn = get_hash_fn(algorithm=algorithm, hash_trim=hash_trim)
    for name, sequence in iter_sequences(fh=fh, fsep=fsep, field_nr=field):
        if name is None:
            continue
        if validate and not all(n in 'ACGT' for n in sequence):
            hexhash = ''
        else:
            data = sequence.strip().upper().encode('ASCII')
            hexhash = hash_fn(data)
            if algprefix:
                hexhash = f'{algorithm}:{hexhash}'
        out.write(f"{name}\t{hexhash}\n")


# -- Run as Script --

def get_argparser():
    parser = argparse.ArgumentParser(
        prog="hash-sequences",
        description="Hash sequences from a multi-FASTA file.",
        epilog="Copyright (C) 2022 Jan Deneweth"
    )
    add_argparser_args(parser=parser)
    return parser


def add_argparser_args(parser: argparse.ArgumentParser):
    parser.add_argument('--algorithm', '-a', default='sha1', type=str, help="The hashing algorithm used")
    parser.add_argument('--hash_trim', '-t', default=None, type=int, help="If specified, the hexadecimal representation of the output hashes are trimmed to the first X characters")
    parser.add_argument('--algprefix', action='store_true', help="Prefix the hash values with the hashing algorithm and a colon")
    parser.add_argument('--validate', default=True, type=argparse.BooleanOptionalAction, help="Whether to validate the sequence as non-ambiguous nucleotides, hash value will be blank when invalid")
    parser.add_argument('--fsep', '-s', default=' ', type=str, help="The character(s) separating fields in the fasta headers")
    parser.add_argument('--field', '-f', default=1, type=int, help="The field in which the sequence name can be found")
    parser.add_argument('--out', '-o', type=argparse.FileType('w'), default=sys.stdout, help="The output destination, defaults to STDOUT")
    parser.add_argument('SEQUENCES', help='The FASTA format input of sequences, defaults to STDIN', type=argparse.FileType('r'), default=sys.stdin)


def main_ns(ns: argparse.Namespace):
    run(
        fh=ns.SEQUENCES,
        out=ns.out,
        algorithm=ns.algorithm,
        hash_trim=ns.hash_trim,
        algprefix=ns.algprefix,
        validate=ns.validate,
        fsep=ns.fsep,
        field=ns.field,
    )


def main(args: list[str] | None = None):
    parser = get_argparser()
    ns = parser.parse_args(args=args)
    main_ns(ns=ns)


if __name__ == '__main__':
    main()


#
#
# END OF FILE
