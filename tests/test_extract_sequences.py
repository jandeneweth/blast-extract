"""
Test extracting sequences
"""

import os
import io
import time
import pytest
import hashlib
import pathlib
import logging
import urllib.request

import blast_extract.extract_sequences
import blast_extract.seqtools


NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TESTDIR_PATH = pathlib.Path(os.path.dirname(__file__))
TESTDATA_PATH = TESTDIR_PATH / "data"


@pytest.fixture(name="download_file")
def get_downloader_func(cache):
    cache_path = cache.mkdir('seq_downloads')

    def download_file(url: str) -> pathlib.Path:
        file_key = f"seq_{hashlib.sha1(url.encode('utf-8')).hexdigest()}"
        file_path = cache_path / file_key
        if not os.path.exists(file_path):
            time.sleep(0.2)  # Bit of rest for servers! (e.g.: NCBI has limit on downloads per second)
            logging.debug(f"Downloading: {url}")
            urllib.request.urlretrieve(url=url, filename=file_path)
        else:
            logging.debug(f"Already got: {url}")
        return file_path

    return download_file


@pytest.fixture(name="retrieve_genome")
def get_genome_retriever(download_file):

    def retrieve_genome(accession: str) -> pathlib.Path:
        url = f"{NCBI_URL}/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
        return download_file(url=url)

    return retrieve_genome


@pytest.fixture(name="achromobacter_dataset")
def achromobacter_dataset(retrieve_genome):
    target_paths = [
        retrieve_genome(accession=accession)
        for accession in [
            "NZ_CP038034.1",
            "NZ_CP043820.1",
            "NZ_CP053986.1",
            "NZ_CP082965.1",
            "NZ_LR134302.1",
        ]
    ]
    references_path = TESTDATA_PATH / "achromobacter_dataset" / "references.fasta"
    return target_paths, references_path


def test_achromobacter_dataset(achromobacter_dataset, tmp_path, capsys):
    target_paths, references_path = achromobacter_dataset
    for target_path in target_paths:
        blast_extract.extract_sequences.main(args=[
            "--dbdir", str(tmp_path / "dbdir"),
            "--references", str(references_path),
            "--pident", "90.0",
            "--pcov", "95.0",
            str(target_path),
        ])
        out, err = capsys.readouterr()
        loci = {
            header.split(' ', 1)[0]
            for header, seq in blast_extract.seqtools.parse_fasta(fh=io.StringIO(out))
        }
        assert len(loci) >= 6
        assert loci.issubset({'rpoB_1', 'nrdA_1', 'nusA_1', 'gltB_1', 'nuoL_1', 'eno_1', 'lepA_1'})


#
#
# END OF FILE
