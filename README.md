# blastMining
============================================================

This program is made to help you mining NCBI BLAST output.

## Requirements

Before able to execute `blastMining`, you need to install the following programs and make sure that
they are executable and available in your `PATH`:

* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

* [TaxonKit - NCBI Taxonomy Toolkit](https://bioinf.shenwei.me/taxonkit/):
	`Please install` this program and download `taxdump`. 
	Follow their [instruction](https://bioinf.shenwei.me/taxonkit/usage/#before-use)  

* [csvtk](https://github.com/shenwei356/csvtk)

* [Python3](https://www.python.org/) 

## Installation
### Option 1
You can easily install this package using [PyPI](https://pypi.org/project/blastMining/)

```bash
$ pip install blastMining
```

### Option 2
Download the latest realese of [blastMining](https://github.com/NuruddinKhoiry/blastMining/releases/download/0.1.0/blastMining-0.1.0.tar.gz) in my Github repository.

Then install it using pip

```bash
$ pip install blastMining-0.1.0.tar.gz
```

## Tutorial
Running blastn
```bash
blastn -query test_data/ASV.fasta -db nt -out test_data/BLASTn.out -outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid" -max_target_seqs 10
```
**Note**: Please strict to the above `blast outfmt`

Next, `mining` your blast result with one of the following methods:

### Method 1. Vote with p. identity cut-off

```bash
$ blastMining vote -i test_data/BLASTn.out -e 1e-03 -n 10 -txl 99,97,95,90,85,80,75 -o Vote_method
```

### Method 2. vote to species level for all

```bash
$ blastMining voteSpecies -i test_data/BLASTn.out -e 1e-03 -n 10 -o VoteSpecies_method
```

### Method 3. LCA 

```bash
$ blastMining lca -i test_data/BLASTn.out -e 1e-03 -n 10 -o lca_method
```

## Command options
```bash
$ blastMining --help

usage: blastMining [-h] [-v] [--debug] {vote,voteSpecies,lca} ...

blastMining 0.1.0

BLAST outfmt 6 only:
("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")

positional arguments:
  {vote,voteSpecies,lca}
    vote                blastMining: voting method with pident cut-off
    voteSpecies         blastMining: vote at species level for all
    lca                 blastMining: lca method

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --debug               debug mode output [False]
```

### Method 1
```bash
$ blastMining vote --help


usage: blastMining vote [-h] -i INPUT [-e EVALUE] [-n TOPN] [-txl TAXA_LEVEL] -o OUTPUT

blastMining: voting method with pident cut-off

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold) [default=1-e3]
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -txl TAXA_LEVEL, --taxa_level TAXA_LEVEL
                        P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
                        [default=99,97,95,90,85,80,75]
  -o OUTPUT, --output OUTPUT
                        output

```

### Method 2
```bash
$ blastMining voteSpecies --help


usage: blastMining voteSpecies [-h] -i INPUT [-e EVALUE] [-n TOPN] -o OUTPUT

blastMining: vote at species level for all

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY: ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold) [default=1-e3]
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -o OUTPUT, --output OUTPUT
                        output
```

### Method 3
```bash
$ blastMining lca --help


usage: blastMining lca [-h] -i INPUT [-e EVALUE] [-n TOPN] -o OUTPUT

blastMining: lca method

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY: ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold) [default=1-e3]
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -o OUTPUT, --output OUTPUT
                        output
```

# Citation
**If you find this package useful**, `please cite` (https://github.com/NuruddinKhoiry/blastMining)