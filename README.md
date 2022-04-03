# Blast-Extract

Extract subsequences from a genome using NCBI BLAST and reference sequences.

Minimal usage:

```bash
extract-sequences --references references.fasta genome.fasta > extracted.fasta
```


## Installation

Requirements:

- [Python 3.10] or later
- NCBI [BLAST+ executables] 2.23.0 or later

The BLAST bin directory must be added to the PATH.

Install the package using pip. This will also create the entry point scripts to use the commands as described below.

```bash
pip install blast-extract
```


## Usage 

All commands can use STDIN and STDOUT, so can be chained using pipes (`|`) on UNIX-like systems.

```bash
cat genome.fasta | extract-sequences --references references.fasta | hash-sequences > seqhashes.tsv
```

### Commands

#### extract-sequences

The `extract-sequences` command extracts nucleotide sequences from an input (multi-)FASTA file. It requires a
multi-FASTA file with references to search for. Its output is a multi-FASTA whose headers start with the identifier
of the corresponding matched reference.

The script will by default normalize and extend the BLAST results so they are in the same strand as the reference and
cover as much of the reference as possible. This is an important improvement over using solely BLAST, as its score
optimization may drop nucleotides at the outer ends of the alignment.

Behind the scenes, a BLAST nucleotides database will be created, by default in the current working directory. You can
specify a different directory for the databases with `--dbdir`. The derived BLAST database will be updated if needed, by
comparing the last modification time of the files. Each reference filepath will generate a uniquely named database, so
can be used in conjunction.

There are two filter options: `--pident` (default 80.0) for minimum percentage identity, and `--pcov` (default 80.0) for
minimum percentage coverage. By nature of using BLAST, less strict settings will increase the likelyhood of
false-negatives: when a subsequence may exist that satisfies the filters, but due to BLAST heuristics was not found. The
underlying blastn call will actually use a slightly lower percentage identity threshold and no coverage requirement.
Both are checked by the python script after post-BLAST corrections.

#### hash-sequences

The `hash-sequences` command is used to get hash values for each sequence in an input multi-FASTA file. You may choose
any of the algorithms available to your python installation, which will be listed if you give an invalid value. It is
recommended to stick to the defaults unless needed for compatibility with other tools or organisations. The output uses
TSV format, with the first column containing the sequence name and the second containing the generated sequence hash.


----

Copyright (C) Jan Deneweth 2022


[Python 3.10]: https://www.python.org/downloads/
[BLAST+ executables]: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
