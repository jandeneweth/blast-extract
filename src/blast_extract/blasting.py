"""
Functionality related to performing a blast.
"""

import os
import sys
import uuid
import subprocess
import typing as t

from blast_extract import blastresult

# We assume BLAST writes to STDOUT with default system encoding?
BLAST_ENCODING = sys.getdefaultencoding()


def run_blast(
        db_basepath: str,
        query: str,
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


def ensure_db(references: str, dbdir: str) -> str:
    """Ensure a database for the given references exists in the database directory, and get the base path of the database files."""
    if not check_db(references=references, dbdir=dbdir):
        make_db(references=references, dbdir=dbdir)
    return get_db_basepath(references=references, dbdir=dbdir)


#
#
# END OF FILE
