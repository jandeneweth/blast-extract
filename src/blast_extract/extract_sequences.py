"""
Extract alleles from a FASTA input file based on reference loci.
"""

import collections
import os
import sys
import argparse
import typing as t

from blast_extract import seqtools
from blast_extract import blasting
from blast_extract import blastresult


def run(
    refspath: str,
    genome: t.TextIO,
    dbdir: str,
    out: t.TextIO,
    fsep: str,
    min_perc_ident: float,
    min_perc_cov: float,
    do_normalize: bool = True,
    do_extend: bool = True,
    filter_best_abbrev: str = 'W',
):
    dbdir = os.path.abspath(os.path.expanduser(os.path.expandvars(dbdir)))
    db_basepath = blasting.ensure_db(refspath=refspath, dbdir=dbdir)
    contig_mapping = seqtools.parse_genome_contigs(genome=genome, fsep=fsep)
    contig_fasta = seqtools.make_contigs_fasta(contigs=contig_mapping.values())
    results = blasting.run_blast(db_basepath=db_basepath, query=contig_fasta, contig_mapping=contig_mapping)
    if do_normalize or do_extend:
        results = (r.normalize() for r in results)
    if do_extend:
        results = (r.extend() for r in results)
    filter_best_func = filter_best_func_map[filter_best_abbrev]
    if filter_best_func:
        results = filter_best_results(results=results, key_func=filter_best_func)
    results = filter_blast_results(results=results, min_perc_ident=min_perc_ident, min_perc_cov=min_perc_cov)
    for result in results:
        out.write(f">{result.blast_header}\n{result.qseq}\n")


def filter_blast_results(results: t.Iterable['blastresult.BlastResult'], min_perc_ident: float, min_perc_cov: float) -> t.Iterable['blastresult.BlastResult']:
    for result in results:
        if result.perc_ident < min_perc_ident:
            continue
        if result.perc_cov < min_perc_cov:
            continue
        yield result


def filter_best_results(
    results: t.Iterable['blastresult.BlastResult'],
    key_func: t.Callable[['blastresult.BlastResult'], float],
) -> list['blastresult.BlastResult']:
    resmap: dict[str, list['blastresult.BlastResult']] = collections.defaultdict(list)
    for result in results:
        resmap[result.sacc].append(result)
    return [max(rs, key=key_func) for rs in resmap.values()]


filter_best_func_map = {
    'N': None,
    'W': (lambda r: (r.perc_ident * r.perc_cov, r.perc_ident)),
    'P': (lambda r: (r.perc_ident, r.perc_cov)),
    'C': (lambda r: (r.perc_cov, r.perc_ident)),
}


# -- Run as Script --

def get_argparser():
    parser = argparse.ArgumentParser(
        prog="extract-sequences",
        description="Extract subsequences from a genome using BLAST.",
        epilog="Copyright (C) 2022 Jan Deneweth"
    )
    add_argparser_args(parser=parser)
    return parser


def add_argparser_args(parser: argparse.ArgumentParser):
    parser.add_argument('--references', '-r', help='The reference alleles FASTA filepath')
    parser.add_argument('--dbdir', '-d', default=os.getcwd(), help="The directory for BLAST databases, defaults to working directory")
    parser.add_argument('--out', '-o', type=argparse.FileType('w'), default=sys.stdout, help="The output destination, defaults to STDOUT")
    parser.add_argument('--fsep', '-s', default=' ', type=str, help="The character(s) separating fields in the fasta headers")
    parser.add_argument('--normalize', default='Y', type=str, choices=['Y', 'N'], help="Normalize results to the 'plus' strand of the reference")
    parser.add_argument('--extend', default='Y', type=str, choices=['Y', 'N'], help="Extend results to cover as much as possible of the reference, implies `normalize`")
    parser.add_argument('--best', default='W', type=str, choices=['N', 'W', 'P', 'C'], help="Whether and how to filter for the best matches. N: No filtering. W: pident*pcov. P: pident. C: coverage.")
    parser.add_argument('--pident', default=80.0, type=float, help="Minimum percent identity for results")
    parser.add_argument('--pcov', default=80.0, type=float, help="Minimum percent coverage for results")
    parser.add_argument('GENOME', help='The FASTA format input assembled genome, defaults to STDIN', type=argparse.FileType('r'), default=sys.stdin)


def main_ns(ns: argparse.Namespace):
    do_extend = True if ns.extend == 'Y' else False
    do_normalize = (True if ns.normalize == 'Y' else False) or do_extend
    run(
        refspath=ns.references,
        genome=ns.GENOME,
        dbdir=ns.dbdir,
        out=ns.out,
        fsep=ns.fsep,
        min_perc_ident=ns.pident,
        min_perc_cov=ns.pcov,
        do_normalize=do_normalize,
        do_extend=do_extend,
        filter_best_abbrev=ns.best,
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
