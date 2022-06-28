# blastMining
============================================================

This program is made to help you mining NCBI BLAST output.

## Requirements

Before able to execute `blastMining`, you need to install the following programs and make sure that
they are executable and available in your `PATH`:

* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

* [TaxonKit - NCBI Taxonomy Toolkit](https://bioinf.shenwei.me/taxonkit/): `Please install` this program and download `taxdump`. Follow their [instruction](https://bioinf.shenwei.me/taxonkit/usage/#before-use)  

* [csvtk](https://github.com/shenwei356/csvtk)

* [Python3](https://www.python.org/) 

**You can install the dependecy programs with conda**

Make sure your conda environment is `up to date` for the sake of the dependency programs.

```bash
$ conda update -n base -c defaults conda

$ conda install -c bioconda blast taxonkit csvtk
```
**Don't forget** to install the required databases for `BLAST` and `TaxonKit`

## Installation

### Option 1

You can easily install this package using [PyPI](https://pypi.org/project/blastMining/)
```bash
$ pip install blastMining
```

### Option 2

Download the latest realese of [blastMining](https://github.com/NuruddinKhoiry/blastMining/releases/download/0.1.1/blastMining-0.1.1.tar.gz) in my Github repository.

Then install it using pip

```bash
$ pip install blastMining-0.1.1.tar.gz
```

## Tutorial
Running blastn
```bash
$ blastn -query test_data/ASV.fasta -db nt -out test_data/BLASTn.out -outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid" -max_target_seqs 10
```
**Note**: Please strict to the above `blast outfmt`

Next, `mining` your blast result with one of the following methods:

### Method 1. Majority vote with percent identity cut-off

The default percent identity cut-off is `99`, `97`, `95`, `90`, `85`, `80`, and `75` for `Species`, `Genus`, `Family`, `Order`, `Class`, `Phylum`, and `Kingdom`, respectively.

```bash
$ blastMining vote -i test_data/BLASTn.out -e 1e-03 -n 10 -txl 99,97,95,90,85,80,75 -sm 'Sample' -o Vote_method
```


### Method 2. Majority vote to species level

```bash
$ blastMining voteSpecies -i test_data/BLASTn.out -e 1e-03 -pi 97 -n 10 -sm 'Sample' -o VoteSpecies_method
```


### Method 3. LCA 

```bash
$ blastMining lca -i test_data/BLASTn.out -e 1e-03 -pi 97 -n 10 -sm 'Sample' -o lca_method
```


### Method 4. besthit 

```bash
$ blastMining besthit -i test_data/BLASTn.out -e 1e-03 -pi 97 -n 10 -sm 'Sample' -o besthit_method
```


## Full_pipeline option

This option allows you to run a full pipeline started from `blastn` -> `blastn_output` ->  `blastMining method` -> `OUTPUT`. 

You can select one of the following combinations:

### BLAST + vote
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m vote -e 1e-03 -txl 99,97,95,90,85,80,75 -n 10 -sm 'Sample' -o vote_pipe
```


### BLAST + voteSpecies
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m voteSpecies -e 1e-03 -pi 97 -n 10 -sm 'Sample' -o voteSpecies_pipe
```


### BLAST + lca
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m lca -e 1e-03 -pi 97 -n 10 -sm 'Sample' -o lca_pipe
```


### BLAST + besthit
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m besthit -e 1e-03 -pi 97 -n 10 -sm 'Sample' -o besthit_pipe
```



## Command options
```bash
$ blastMining --help


usage: blastMining [-h] [-v] [--debug] {vote,voteSpecies,lca,besthit,full_pipeline} ...

blastMining 0.1.1

BLAST outfmt 6 only:
("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")

positional arguments:
  {vote,voteSpecies,lca,besthit,full_pipeline}
    vote                blastMining: voting method with pident cut-off
    voteSpecies         blastMining: vote at species level for all
    lca                 blastMining: lca method
    besthit             blastMining: besthit method
    full_pipeline       blastMining: Running BLAST + mining the output

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --debug               debug mode output [False]
```



### Method 1
```bash
$ blastMining vote --help


usage: blastMining vote [-h] -i INPUT [-e EVALUE] [-txl TAXA_LEVEL] [-n TOPN] [-sm SAMPLE_NAME] -o OUTPUT

blastMining: voting method with pident cut-off

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY: ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold) [default=1-e3]
  -txl TAXA_LEVEL, --taxa_level TAXA_LEVEL
                        P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species [default=99,97,95,90,85,80,75]
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -o OUTPUT, --output OUTPUT
                        output
```



### Method 2
```bash
$ blastMining voteSpecies --help


usage: blastMining voteSpecies [-h] -i INPUT [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] -o OUTPUT

blastMining: vote at species level for all

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY: ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold) [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity (Ignore hits if their p. identities are below this threshold) [default=97]
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -o OUTPUT, --output OUTPUT
                        output
```



### Method 3
```bash
$ blastMining lca --help


usage: blastMining lca [-h] -i INPUT [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] -o OUTPUT

blastMining: lca method

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY: ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold) [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity (Ignore hits if their p. identities are below this threshold) [default=97]
  -n TOPN, --topN TOPN  Top N hits used for LCA calculation [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -o OUTPUT, --output OUTPUT
                        output
```



### Method 4
```bash
$ blastMining besthit --help


usage: blastMining besthit [-h] -i INPUT [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] -o OUTPUT

blastMining: besthit method

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY: ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold) [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity (Ignore hits if their p. identities are below this threshold) [default=97]
  -n TOPN, --topN TOPN  Top N hits used for sorting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -o OUTPUT, --output OUTPUT
                        output
```



### Full pipeline
```bash
$ full_pipeline --help


usage: blastMining full_pipeline [-h] -i INPUT -bp BLAST_PARAM [-m MINING] [-e EVALUE] [-pi PIDENT] [-txl TAXA_LEVEL] [-n TOPN] [-sm SAMPLE_NAME] -o OUTPUT

blastMining: Running BLAST + mining the output

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input FASTA
  -bp BLAST_PARAM, --blast_param BLAST_PARAM
                        BLAST parameters: Note: "-outfmt" has been defined by the package, you don't need to add it (-outfmt="6 qseqid sseqid pident length mismatch
                        gapopen evalue bitscore staxid") [default="-db nt -num_threads 1 -max_target_seqs 10"]
  -m MINING, --mining MINING
                        blast mining method Available methods={'vote','voteSpecies','lca','besthit'} [default='vote']
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue (Ignore hits if their evalues are above this threshold ) [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity (Ignore hits if their p. identities are below this threshold) [default=97] **Not compatible** with "vote method"
  -txl TAXA_LEVEL, --taxa_level TAXA_LEVEL
                        P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species [default=99,97,95,90,85,80,75] **Required** for "vote method"
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -o OUTPUT, --output OUTPUT
                        output
```


# Citation
**If you find this package useful**, `please cite` (https://github.com/NuruddinKhoiry/blastMining)
