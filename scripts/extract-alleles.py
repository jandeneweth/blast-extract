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

import Bio.Seq
import Bio.SeqIO


# We assume BLAST writes to STDOUT with default system encoding?
BLAST_ENCODING = sys.getdefaultencoding()
GAPCHAR = '-'  # Assume dash used for length-1 gap (in IUPAC actually specified as undefined length)
# Note: 'qstrand' doesn't exist in BLAST... we'll have to derive it from start and end positions.
BLAST_OUTFIELDS = ['qacc', 'sacc', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'qseq', 'sseq', 'slen']

BlastResult = collections.namedtuple('BlastResult', BLAST_OUTFIELDS)


def run(
    references: str,
    genome: typing.TextIO,
    dbdir: str,
    out: typing.TextIO,
    fsep: str,
    min_perc_ident: float,
    min_perc_cov: float,
):
    dbdir = os.path.abspath(os.path.expanduser(os.path.expandvars(dbdir)))
    if not check_db(references=references, dbdir=dbdir):
        make_db(references=references, dbdir=dbdir)
    db_basepath = get_db_basepath(references=references, dbdir=dbdir)
    contig_mapping = parse_genome_contigs(genome=genome, fsep=fsep)
    contig_fasta = make_contigs_fasta(contig_mapping=contig_mapping)
    results = list(run_blast(query=contig_fasta, db_basepath=db_basepath))
    results = iter_norm_extend_seqs(contig_mapping=contig_mapping, results=results)
    results = filter_blast_results(results=results, min_perc_ident=min_perc_ident, min_perc_cov=min_perc_cov)
    for result in results:
        out.write(f">{result.sacc}\n{result.qseq}\n")


def iter_norm_extend_seqs(contig_mapping: dict[str, str], results: typing.Iterable['BlastResult']) -> typing.Iterable:
    """Normalize strand and extend query sequences."""
    for result in results:
        assert result.sstrand in ('minus', 'plus'), f"Unexpected sstrand: '{result.sstrand}'"
        if result.sstrand == 'minus':  # On reverse complement => normalize
            result = _revcomp_result(result=result)
        qsrcseq = contig_mapping[result.qacc]
        if result.sstart > 1:  # Missing start
            result = _extend_result_qry_start(
                result=result,
                positions=result.sstart-1,
                qsrcseq=qsrcseq
            )
        if result.send < result.slen:  # Missing end
            result = _extend_result_qry_end(
                result=result,
                positions=result.slen-result.send,
                qsrcseq=qsrcseq
            )
        yield result


def _revcomp_result(result: 'BlastResult') -> 'BlastResult':
    new_strand = 'plus' if result.sstrand == 'minus' else 'minus'
    new_qseq = Bio.Seq.reverse_complement(result.qseq)
    new_sseq = Bio.Seq.reverse_complement(result.sseq)
    return BlastResult(
        qacc=result.qacc,
        sacc=result.sacc,
        qstart=result.qend,  # Swapped
        qend=result.qstart,  # Swapped
        sstart=result.send,  # Swapped
        send=result.sstart,  # Swapped
        sstrand=new_strand,  # Switched
        qseq=new_qseq,  # Reverse complemented
        sseq=new_sseq,  # Reverse complemented
        slen=result.slen
    )


def _extend_result_qry_start(result: 'BlastResult', positions: int, qsrcseq: str):
    assert result.sstrand == 'plus', "Only use normalized references!"
    # Reverse complement the query source sequence if needed
    qstrand = 'plus' if result.qend >= result.qstart else 'minus'
    if qstrand == 'minus':
        qsrcseq = Bio.Seq.reverse_complement(qsrcseq)
    # Determine maximal extendability and extension sequence
    if qstrand == 'plus':
        if positions > (result.qstart - 1):
            positions = result.qstart - 1
        # Account for python base-0 coord system, substract one from start
        qsrcend = result.qstart-1
        qext = qsrcseq[qsrcend-positions:qsrcend]
        sext = '-' * positions
        new_qstart = result.qstart - positions
    elif qstrand == 'minus':
        if positions > (len(qsrcseq) - result.qstart):
            positions = len(qsrcseq) - result.qstart
        # Ex.: Start pos 48 in len 50 seq is pos 3 in reverse complemented form (base-1 coord system), but substract one for python base-0 coord system => pos 2
        qsrcend = len(qsrcseq) - result.qstart
        qext = qsrcseq[qsrcend-positions:qsrcend]
        sext = '-' * positions
        new_qstart = result.qstart + positions
    else:
        raise AssertionError("Unhandled strand")
    # Extend query and ref sequences, the latter one with gaps (dashes)
    new_qseq = qext + result.qseq
    new_sseq = sext + result.sseq
    # Return new result
    return BlastResult(
        qacc=result.qacc,
        sacc=result.sacc,
        qstart=new_qstart,  # Moved 'up'
        qend=result.qend,
        sstart=result.sstart - positions,  # Moved 'up'
        send=result.send,
        sstrand=result.sstrand,
        qseq=new_qseq,  # Extended
        sseq=new_sseq,  # Extended
        slen=result.slen
    )


def _extend_result_qry_end(result: 'BlastResult', positions: int, qsrcseq: str):
    assert result.sstrand == 'plus', "Only use normalized references!"
    # Reverse complement the query source sequence if needed
    qstrand = 'plus' if result.qend >= result.qstart else 'minus'
    if qstrand == 'minus':
        qsrcseq = Bio.Seq.reverse_complement(qsrcseq)
    # Determine maximal extendability and extension sequence
    if qstrand == 'plus':
        if positions > (len(qsrcseq) - result.qend):
            positions = len(qsrcseq) - result.qend
        # Account for python base-0 coord system, don't need to +1 for skipping a character
        qsrcstart = result.qend
        qext = qsrcseq[qsrcstart:qsrcstart+positions]
        sext = '-' * positions
        new_qend = result.qend + positions
    elif qstrand == 'minus':
        if positions > (result.qstart - 1):
            positions = result.qstart - 1
        # Ex.: End pos 5 in len 50 seq is pos 45 in reverse complemented form (base-1 coord system), still need to +1 for exclusion in base-0 coord system
        qsrcstart = len(qsrcseq) - result.qend + 1
        qext = qsrcseq[qsrcstart:qsrcstart+positions]
        sext = '-' * positions
        new_qend = result.qend - positions
    else:
        raise AssertionError("Unhandled strand")
    # Extend query and ref sequences, the latter one with gaps (dashes)
    new_qseq = qext + result.qseq
    new_sseq = sext + result.sseq
    # Return new result
    return BlastResult(
        qacc=result.qacc,
        sacc=result.sacc,
        qstart=result.qstart,
        qend=new_qend,  # Moved 'down'
        sstart=result.sstart,
        send=result.send + positions,  # Moved 'down'
        sstrand=result.sstrand,
        qseq=new_qseq,  # Extended
        sseq=new_sseq,  # Extended
        slen=result.slen
    )


def filter_blast_results(results: typing.Iterable['BlastResult'], min_perc_ident: float, min_perc_cov: float) -> typing.Iterable['BlastResult']:
    for result in results:
        perc_ident, perc_cov = _calc_perc_ident_cov(qry=result.qseq, subj=result.sseq)
        if perc_ident < min_perc_ident:
            continue
        if perc_cov < min_perc_cov:
            continue
        yield result


def _calc_perc_ident_cov(qry: str, subj: str):
    # Assuming nucleotide sequences, 1 for exact equals, 0 otherwise (no intermediate for IUPAC).
    # Assuming both sequences are equal length, strict zip will raise error otherwise.
    # Assuming qry and subj never both have a gap character at the same position.
    ident_length = len(qry)
    cov_length = 0
    ident_score = 0
    cov_score = 0
    for qry_c, subj_c in zip(qry, subj, strict=True):
        ident_score += qry_c == subj_c
        cov_score += qry_c != GAPCHAR and subj_c != GAPCHAR  # Note: Only count coverage where subject does not have gap (=> insert may not compensate deletion!)
        cov_length += subj_c != GAPCHAR
    perc_ident = ident_score / ident_length * 100.0
    perc_cov = cov_score / cov_length * 100.0
    return perc_ident, perc_cov


def make_contigs_fasta(contig_mapping: dict[str, str]) -> str:
    result = ''
    for name, seq in contig_mapping.items():
        result += f">{name}\n{seq}\n"
    return result


def parse_genome_contigs(genome: typing.TextIO, fsep) -> dict[str, str]:
    """Parse a genome fasta file into a mapping of contig names and sequences."""
    contig_mapping = dict()
    for seqrecord in Bio.SeqIO.parse(genome, 'fasta'):
        name = seqrecord.description.split(fsep, 1)[0]  # First value in header assumed to be contig name
        seq = str(seqrecord.seq).lower()
        contig_mapping[name] = seq
    return contig_mapping


def run_blast(
        query: str,
        db_basepath: str,
        perc_ident_guideline: float | None = None,
) -> typing.Iterable['BlastResult']:
    args = ['blastn', '-db', db_basepath, '-outfmt', f'6 {" ".join(BLAST_OUTFIELDS)}', ]
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
    results = []
    for line in result.stdout.splitlines():
        fields = line.split('\t')
        raw = BlastResult(*fields)
        result = BlastResult(
            qacc=raw.qacc,
            sacc=raw.sacc,
            qstart=int(raw.qstart),
            qend=int(raw.qend),
            sstart=int(raw.sstart),
            send=int(raw.send),
            sstrand=raw.sstrand,
            qseq=raw.qseq,
            sseq=raw.sseq,
            slen=int(raw.slen),
        )
        results.append(result)
    return results


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


def main(args: list[str]):
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
    main(args=sys.argv[1:])

#
#
# END OF FILE
