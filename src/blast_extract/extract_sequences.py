#!/usr/bin/env python3
"""
Extract alleles from a FASTA input file based on reference loci.
"""

import os
import sys
import uuid
import argparse
import subprocess
import typing as t

import Bio.Seq
import Bio.SeqIO

from blast_extract import blastresult


# We assume BLAST writes to STDOUT with default system encoding?
BLAST_ENCODING = sys.getdefaultencoding()


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
    if not check_db(references=references, dbdir=dbdir):
        make_db(references=references, dbdir=dbdir)
    db_basepath = get_db_basepath(references=references, dbdir=dbdir)
    contig_mapping = parse_genome_contigs(genome=genome, fsep=fsep)
    contig_fasta = make_contigs_fasta(contigs=contig_mapping.values())
    results = list(run_blast(query=contig_fasta, db_basepath=db_basepath, contig_mapping=contig_mapping))
    results = [r.normalize_and_extend() for r in results]
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


def run_blast(
        query: str,
        db_basepath: str,
        contig_mapping: dict[str, 'blastresult.Contig'],
        perc_ident_guideline: float | None = None,
) -> t.Iterable['blastresult.BlastResult']:
    args = ['blastn', '-db', db_basepath, '-outfmt', blastresult.BlastResult.REQUIRED_OUTFMT, ]
    if perc_ident_guideline:
        # Note: set 2% lower, since there may be extension before further filtering, which may influence the value
        args.extend(['--perc_ident', max(0.0, perc_ident_guideline - 2.0)])
    result = subprocess.run(
        args=args,
        input=query,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding=BLAST_ENCODING,
    )
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(f"Error running BLAST: {result.stderr or '<no stderr>'} (returncode {result.returncode})")
    return blastresult.BlastResult.from_output(output=result.stdout, contig_mapping=contig_mapping)


def make_db(references: str, dbdir: str):
    """Create the BLAST database for the references."""
    db_basepath = get_db_basepath(references=references, dbdir=dbdir)
    db_pardir = os.path.dirname(db_basepath)
    os.makedirs(db_pardir, exist_ok=True)
    args = ['makeblastdb', '-dbtype', 'nucl', '-out', db_basepath, '-in', references]
    result = subprocess.run(
        args=args,
        stderr=subprocess.PIPE,
        encoding=BLAST_ENCODING,
    )
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(f"Error building BLAST database: {result.stderr or '<no stderr>'} (returncode {result.returncode})")


def check_db(references: str, dbdir: str) -> bool:
    """Whether the BLAST database for the references exists and is up-to-date."""
    db_basepath = get_db_basepath(references=references, dbdir=dbdir)
    testpath = f"{db_basepath}.nhr"  # The 'header' file
    try:
        db_mtime = os.path.getmtime(testpath)
    except OSError:
        return False  # Does not exist
    else:
        try:
            refs_mtime = os.path.getmtime(references)
        except OSError as e:
            raise RuntimeError(f"Invalid references filepath '{references}': {e}")
        if refs_mtime > db_mtime:
            return False  # DB is outdated
    return True


def get_db_basepath(references: str, dbdir: str) -> str:
    """Get a unique, reproducible, safe database path for a given filepath."""
    unique = str(uuid.uuid5(namespace=uuid.NAMESPACE_URL, name=references))
    basename = f"db_{unique}"  # Prefixed to ensure it doesn't start with a digit
    return os.path.join(dbdir, basename)


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
