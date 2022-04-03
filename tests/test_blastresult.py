"""
BlastResult tests
"""

from blast_extract import seqtools
from blast_extract import blastresult


def test_normalize():
    result = blastresult.BlastResult(
        qacc='testqry',
        sacc='testsub',
        qstart=1,
        qend=3,
        sstart=3,
        send=1,
        qseq="AAA",
        sseq="AAA",
        slen=3,
        qcontig=seqtools.Sequence(name='testqry', fwdseq='AAA'),
    )
    result = result.normalize()
    assert result.qseq == "TTT"
    assert result.sseq == "TTT"
    assert result.qstart == 3
    assert result.qend == 1
    assert result.sstart == 1
    assert result.send == 3


def test_extend():
    result = _make_result(sstart=1+2, send=_TESTSUB_LEN-2)  # Missing 2 at start and end
    result = result.extend()
    assert len(result.sseq) == _TESTSUB_LEN
    assert result.qseq == _TESTSUB_REAL  # Equal because no differences between qry and subj where aligned
    assert result.sstart == 1
    assert result.send == _TESTSUB_LEN


_TESTQRY = "CCCGGGATGCCCGGGCCCGGGCCCTAGCCCGGG"
_TESTSUB = "      ATGCCCGGGCCCGGGCCCTAG      "
_TESTSUB_REAL = _TESTSUB.replace(' ', '')
_TESTSUB_LEN = len(_TESTSUB_REAL)
_ADD_START = 6


def _make_result(sstart, send):
    sseq = _TESTSUB[sstart+_ADD_START-1:send+_ADD_START]
    qseq = _TESTQRY[sstart+_ADD_START-1:send+_ADD_START]
    qstart = sstart+_ADD_START
    qend = send+_ADD_START
    qcontig = seqtools.Sequence(name='testqry', fwdseq=_TESTQRY)
    return blastresult.BlastResult(
        qacc=qcontig.name,
        sacc='testsub',
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=send,
        qseq=qseq,
        sseq=sseq,
        slen=_TESTSUB_LEN,
        qcontig=qcontig
    )


#
#
# END OF FILE
