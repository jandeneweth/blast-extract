"""
Blast extract CLI
"""

import argparse

from blast_extract import hash_sequences
from blast_extract import extract_sequences


def get_argparser():
    parser = argparse.ArgumentParser(
        prog="blast-extract",
        description="Blast extraction package CLI.",
        epilog="Copyright (C) 2022 Jan Deneweth"
    )
    subparsers = parser.add_subparsers(title='Subcommands', metavar='SUBCOMMAND', required=True)
    hash_parser = subparsers.add_parser(name='hash-sequences', help="Hash sequences from a multi-FASTA file.")
    hash_parser.set_defaults(subcommand='hash-sequences')
    hash_sequences.add_argparser_args(parser=hash_parser)
    extract_parser = subparsers.add_parser(name='extract-sequences', help="Extract subsequences from a genome using BLAST.")
    extract_parser.set_defaults(subcommand='extract-sequences')
    extract_sequences.add_argparser_args(parser=extract_parser)
    return parser


def main(args: list[str] | None = None):
    parser = get_argparser()
    ns = parser.parse_args(args=args)
    if ns.subcommand == 'hash-sequences':
        hash_sequences.main_ns(ns=ns)
    elif ns.subcommand == 'extract-sequences':
        extract_sequences.main_ns(ns=ns)
    else:
        raise AssertionError(f"Unhandled subcommand: '{ns.subcommand}'")


if __name__ == '__main__':
    main()


#
#
# END OF FILE
