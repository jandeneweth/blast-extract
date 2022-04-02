"""
Blast result related objects.
"""

import typing as t
import collections

from blast_extract import seqtools

GAPCHAR = '-'  # Assume dash used for length-1 gap (in IUPAC actually specified as undefined length)
_BLAST_OUTFIELDS = ['qacc', 'sacc', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'slen']
_RawBlastResult = collections.namedtuple('BlastResult', _BLAST_OUTFIELDS)


class BlastResult:

    REQUIRED_OUTFMT = f'6 {" ".join(_BLAST_OUTFIELDS)}'

    def __init__(
            self,
            qacc: str,
            sacc: str,
            qstart: int,
            qend: int,
            sstart: int,
            send: int,
            qseq: str,
            sseq: str,
            slen: int,
            qcontig: 'seqtools.Contig',
    ):
        self.qacc = qacc
        self.sacc = sacc
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.qseq = qseq.upper()
        self.sseq = sseq.upper()
        self.slen = slen
        self.qcontig = qcontig
        self._perc_ident = None
        self._perc_cov = None

    @property
    def qstrand(self):
        return 'plus' if self.qend >= self.qstart else 'minus'

    @property
    def sstrand(self):
        return 'plus' if self.send >= self.sstart else 'minus'

    @property
    def perc_ident(self):
        """Percentage identity"""
        if self._perc_ident is None:
            self._perc_ident, self._perc_cov = self._calc_perc_ident_cov()
        return self._perc_ident

    @property
    def perc_cov(self):
        """Percentage coverage"""
        if self._perc_cov is None:
            self._perc_ident, self._perc_cov = self._calc_perc_ident_cov()
        return self._perc_cov

    @classmethod
    def from_output(cls, output: str, contig_mapping: dict[str, 'seqtools.Contig']) -> t.Generator['BlastResult', None, None]:
        """Create results from blastn output."""
        for line in output.splitlines():
            fields = line.split('\t')
            raw = _RawBlastResult(*fields)
            result = cls(
                qacc=raw.qacc,
                sacc=raw.sacc,
                qstart=int(raw.qstart),
                qend=int(raw.qend),
                sstart=int(raw.sstart),
                send=int(raw.send),
                qseq=raw.qseq,
                sseq=raw.sseq,
                slen=int(raw.slen),
                qcontig=contig_mapping[raw.qacc]
            )
            yield result

    def normalize_and_extend(self) -> 'BlastResult':
        """Normalize strand and extend the sequences to cover as much of the full subject sequence as possible."""
        result = self
        if self.sstrand == 'minus':  # On reverse complement => normalize
            result = self._revcomp_result(result=result)
        if result.sstart > 1:  # Missing start
            result = self._extend_result_qry_start(
                result=result,
                positions=result.sstart-1,
            )
        if result.send < result.slen:  # Missing end
            result = self._extend_result_qry_end(
                result=result,
                positions=result.slen-result.send,
            )
        return result

    def _calc_perc_ident_cov(self):
        # Assuming nucleotide sequences, 1 for exact equals, 0 otherwise (no intermediate for IUPAC).
        # Assuming both sequences are equal length, strict zip will raise error otherwise.
        # Assuming qry and subj never both have a gap character at the same position.
        ident_length = len(self.qseq)
        cov_length = 0
        ident_score = 0
        cov_score = 0
        for qry_c, subj_c in zip(self.qseq, self.sseq, strict=True):
            ident_score += qry_c == subj_c
            cov_score += qry_c != GAPCHAR and subj_c != GAPCHAR  # Note: Only count coverage where subject does not have gap (=> insert may not compensate deletion!)
            cov_length += subj_c != GAPCHAR
        perc_ident = ident_score / ident_length * 100.0
        perc_cov = cov_score / cov_length * 100.0
        return perc_ident, perc_cov

    @staticmethod
    def _revcomp_result(result: 'BlastResult') -> 'BlastResult':
        new_qseq = seqtools.reverse_complement(result.qseq)
        new_sseq = seqtools.reverse_complement(result.sseq)
        return BlastResult(
            qacc=result.qacc,
            sacc=result.sacc,
            qstart=result.qend,  # Swapped
            qend=result.qstart,  # Swapped
            sstart=result.send,  # Swapped
            send=result.sstart,  # Swapped
            qseq=new_qseq,  # Reverse complemented
            sseq=new_sseq,  # Reverse complemented
            slen=result.slen,
            qcontig=result.qcontig,
        )

    @staticmethod
    def _extend_result_qry_start(result: 'BlastResult', positions: int):
        assert result.sstrand == 'plus', "Only use normalized references!"
        # Determine maximal extendability and extension sequence
        if result.qstrand == 'plus':
            if positions > (result.qstart - 1):
                positions = result.qstart - 1
            # Account for python base-0 coord system, substract one from start
            qsrcend = result.qstart-1
            qext = result.qcontig.fwdseq[qsrcend-positions:qsrcend]
            sext = '-' * positions
            new_qstart = result.qstart - positions
        elif result.qstrand == 'minus':
            if positions > (result.qcontig.length - result.qstart):
                positions = result.qcontig.length - result.qstart
            # Ex.: Start pos 48 in len 50 seq is pos 3 in reverse complemented form (base-1 coord system), but substract one for python base-0 coord system => pos 2
            qsrcend = result.qcontig.length - result.qstart
            qext = result.qcontig.revseq[qsrcend-positions:qsrcend]
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
            qseq=new_qseq,  # Extended
            sseq=new_sseq,  # Extended
            slen=result.slen,
            qcontig=result.qcontig,
        )

    @staticmethod
    def _extend_result_qry_end(result: 'BlastResult', positions: int):
        assert result.sstrand == 'plus', "Only use normalized references!"
        # Determine maximal extendability and extension sequence
        if result.qstrand == 'plus':
            if positions > (result.qcontig.length - result.qend):
                positions = result.qcontig.length - result.qend
            # Account for python base-0 coord system, don't need to +1 for skipping a character
            qsrcstart = result.qend
            qext = result.qcontig.fwdseq[qsrcstart:qsrcstart+positions]
            sext = '-' * positions
            new_qend = result.qend + positions
        elif result.qstrand == 'minus':
            if positions > (result.qstart - 1):
                positions = result.qstart - 1
            # Ex.: End pos 5 in len 50 seq is pos 45 in reverse complemented form (base-1 coord system), still need to +1 for exclusion in base-0 coord system
            qsrcstart = result.qcontig.length - result.qend + 1
            qext = result.qcontig.revseq[qsrcstart:qsrcstart+positions]
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
            qseq=new_qseq,  # Extended
            sseq=new_sseq,  # Extended
            slen=result.slen,
            qcontig=result.qcontig,
        )


#
#
# END OF FILE
