#!/usr/bin/env python3
"""
Extract alleles from a FASTA input file based on reference loci.
"""

import os
import sys
import uuid
import typing
import argparse
import subprocess
import collections


BLAST_ENCODING = sys.getdefaultencoding()  # We assume BLAST writes to STDOUT with default system encoding?
BLAST_OUTFIELDS = ['sacc', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'qseq', 'sseq', ]
BlastResult = collections.namedtuple('BlastResult', BLAST_OUTFIELDS)


def run(
    references: str,
    genome: typing.TextIO,
    dbdir: str,
    out: typing.TextIO,
):
    dbdir = os.path.abspath(os.path.expanduser(os.path.expandvars(dbdir)))
    if not check_db(references=references, dbdir=dbdir):
        make_db(references=references, dbdir=dbdir)
    db_basepath = get_db_basepath(references=references, dbdir=dbdir)
    results = list(run_blast(query=genome, db_basepath=db_basepath))
    # TODO: filter and refine results
    # TODO: write to output FASTA file


def run_blast(query: typing.TextIO, db_basepath: str) -> typing.Iterable['BlastResult']:
    result = subprocess.run(
        args=['blastn', '-db', db_basepath, '-outfmt', f'6 {" ".join(BLAST_OUTFIELDS)}', ],
        stdin=query,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding=BLAST_ENCODING,
    )
    if result.returncode != 0 or result.stderr:
        raise RuntimeError(f"Error running BLAST: {result.stderr or '<no stderr>'} (returncode {result.returncode})")
    results = []
    for line in result.stdout.splitlines():
        fields = line.split('\t')
        results.append(BlastResult(*fields))
    return results


def make_db(references: str, dbdir: str):
    """Create the BLAST database for the references."""
    db_basepath = get_db_basepath(references=references, dbdir=dbdir)
    db_pardir = os.path.dirname(db_basepath)
    os.makedirs(db_pardir, exist_ok=True)
    args = ['makeblastdb', '-dbtype', 'nucl', '-out', db_basepath, '-in', references]
    print(args)
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
    parser.add_argument('genome', help='The FASTA format input assembled genome, defaults to STDIN', type=argparse.FileType('r'), default=sys.stdin)
    return parser


def main(args: list[str]):
    print(args)
    parser = get_argparser()
    ns = parser.parse_args(args=args)
    run(
        references=ns.references,
        genome=ns.genome,
        dbdir=ns.dbdir,
        out=ns.out
    )


if __name__ == '__main__':
    main(args=sys.argv[1:])

#
#
# END OF FILE
