#!/usr/bin/env python3
"""
Extract alleles from a FASTA input file based on reference loci.
"""

import os
import sys
import argparse
import typing as t

import Bio.Seq
import Bio.SeqIO

from blast_extract import blasting
from blast_extract import blastresult


def run(
    references: str,
    genome: t.TextIO,
    dbdir: str,
    out: t.TextIO,
    fsep: str,
    min_perc_ident: float,
    min_perc_cov: float,
):
    dbdir = os.path.abspath(os.path.expanduser(os.path.expandvars(dbdir)))
    db_basepath = blasting.ensure_db(references=references, dbdir=dbdir)
    contig_mapping = parse_genome_contigs(genome=genome, fsep=fsep)
    contig_fasta = make_contigs_fasta(contigs=contig_mapping.values())
    results = blasting.run_blast(db_basepath=db_basepath, query=contig_fasta, contig_mapping=contig_mapping)
    results = (r.normalize_and_extend() for r in results)
    results = filter_blast_results(results=results, min_perc_ident=min_perc_ident, min_perc_cov=min_perc_cov)
    for result in results:
        out.write(f">{result.sacc}\n{result.qseq}\n")


def filter_blast_results(results: t.Iterable['blastresult.BlastResult'], min_perc_ident: float, min_perc_cov: float) -> t.Iterable['blastresult.BlastResult']:
    for result in results:
        if result.perc_ident < min_perc_ident:
            continue
        if result.perc_cov < min_perc_cov:
            continue
        yield result


def make_contigs_fasta(contigs: t.Iterable['blastresult.Contig']) -> str:
    result = ''
    for contig in contigs:
        result += f">{contig.name}\n{contig.fwdseq}\n"
    return result


def parse_genome_contigs(genome: t.TextIO, fsep) -> dict[str, 'blastresult.Contig']:
    """Parse a genome fasta file into a mapping of contig names and sequences."""
    contig_mapping = dict()
    for seqrecord in Bio.SeqIO.parse(genome, 'fasta'):
        name = seqrecord.description.split(fsep, 1)[0]  # First value in header assumed to be contig name
        seq = str(seqrecord.seq).lower()
        contig_mapping[name] = blastresult.Contig(name=name, fwdseq=seq)
    return contig_mapping


# -- Run as Script --

def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--references', '-r', help='The reference alleles FASTA file')
    parser.add_argument('--dbdir', '-d', default=os.getcwd(), help="The directory for BLAST databases, defaults to working directory")
    parser.add_argument('--out', '-o', type=argparse.FileType('w'), default=sys.stdout, help="The output destination, defaults to STDOUT")
    parser.add_argument('--fsep', '-s', default=' ', type=str, help="The character(s) separating fields in the fasta headers")
    parser.add_argument('--pident', default=80.0, type=float, help="Minimum percent identity for results")
    parser.add_argument('--pcov', default=80.0, type=float, help="Minimum percent coverage for results")
    parser.add_argument('genome', help='The FASTA format input assembled genome, defaults to STDIN', type=argparse.FileType('r'), default=sys.stdin)
    return parser


def main(args: list[str] | None = None):
    parser = get_argparser()
    ns = parser.parse_args(args=args)
    run(
        references=ns.references,
        genome=ns.genome,
        dbdir=ns.dbdir,
        out=ns.out,
        fsep=ns.fsep,
        min_perc_ident=ns.pident,
        min_perc_cov=ns.pcov
    )


if __name__ == '__main__':
    main()

#
#
# END OF FILE
