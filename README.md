# blastMining

### Mining BLAST OUTPUT

[![PyPI Version](https://img.shields.io/pypi/v/blastMining?style=flat-square)](https://pypi.org/project/blastMining)
[![Conda Version](https://anaconda.org/bioconda/blastmining/badges/version.svg)](https://anaconda.org/bioconda/blastmining)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.6823244.svg?style=flat-square)](https://zenodo.org/record/6823244)

---

`blastMining` is a tool used for mining NCBI BLAST output from a single or multiple sequences,
including but not limited to ASV/OTU from amplicon sequencing, 
contigs/scaffolds from shotgun metagenomics, etc.

`blastMining` is written in Python (tested with v3.6+). It is
available on the [Python Package Index](https://pypi.org/project/blastMining/)

## Requirements

Before able to execute `blastMining`, you need to install the following programs and make sure that
they are executable and available in your `PATH`:

* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

* [TaxonKit - NCBI Taxonomy Toolkit](https://bioinf.shenwei.me/taxonkit/): `Please install` this program and download `taxdump`. Follow their [instruction](https://bioinf.shenwei.me/taxonkit/usage/#before-use)  

* [csvtk](https://github.com/shenwei356/csvtk)

* [Python3](https://www.python.org/)

* [krona](https://github.com/marbl/Krona/wiki) 


## Installation

### Option 1. Install via [conda](https://anaconda.org/bioconda/blastmining) 
This option will automatically install the dependecy programs. So, you don't need to install them manually. 

```bash
$ conda install -c bioconda blastmining
```
In the case your `blast version` is lower than `2.12.0`, you can update it using `conda`.
```bash
$ conda install -c bioconda blast=2.12.0
```

### Option 2. Install via [PyPI](https://pypi.org/project/blastMining/)

```bash
$ pip install blastMining
```

### Option 3. Install manually

Download the latest realese of [blastMining](https://github.com/NuruddinKhoiry/blastMining/releases/download/1.0.0/blastMining-1.0.0.tar.gz) in my Github repository.

Then install it using pip

```bash
$ pip install blastMining-1.0.0.tar.gz
```

### Installation Notes

If you install `blastMining` using **option 2** or **option 3**, you need to install the dependency programs.

**You can install the dependecy programs with conda**

Make sure your conda environment is `up to date` for the sake of the dependency programs.

```bash
$ conda update -n base conda

$ conda install -c bioconda taxonkit csvtk krona blast=2.12.0
```

### Before use

**Don't forget** to install the required databases for `BLAST` and `TaxonKit`


## Tutorial
Running blastn
```bash
$ blastn -query test_data/ASV.fasta -db nt -out test_data/BLASTn.out -outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid" -max_target_seqs 10
```
**Note**: Please strict to the above `blast outfmt`

Next, `mining` your blast result with **one of the following methods**:

### * Method A. Majority vote with percent identity cut-off

The `vote` algorithm is as follow:

![vote method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/vote_method.png?raw=true)

The default percent identity cut-off is `99`, `97`, `95`, `90`, `85`, `80`, and `75` for `Species`, `Genus`, `Family`, `Order`, `Class`, `Phylum`, and `Kingdom`, respectively.

```bash
$ blastMining vote -i test_data/BLASTn.out -e 1e-03 -n 10 -txl 99,97,95,90,85,80,75 -sm 'Sample' -kp -o Vote_method
```


### * Method B. Majority vote at species level

The `voteSpecies` algorithm is as follow:

![vote method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/voteSpecies_method.png?raw=true)


```bash
$ blastMining voteSpecies -i test_data/BLASTn.out -e 1e-03 -pi 97 -n 10 -sm 'Sample' -kp -o VoteSpecies_method
```


### * Method C. LCA

The `lca` algorithm in is as follow:

![lca method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/lca_method.png?raw=true) 

The `lca` algorithm used in `blastMining` is from [TaxonKit](https://bioinf.shenwei.me/taxonkit/usage/#lca).

```bash
$ blastMining lca -i test_data/BLASTn.out -e 1e-03 -pi 97 -n 10 -sm 'Sample' -kp -o lca_method
```


### * Method D. besthit

The `besthit` algorithm is as follow:

![besthit method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/besthit_method.png?raw=true) 

```bash
$ blastMining besthit -i test_data/BLASTn.out -e 1e-03 -pi 97 -n 10 -sm 'Sample' -kp -o besthit_method
```

----

## Full_pipeline option

This option allows you to run a full pipeline started from `blastn` -> `blastn_output` ->  `blastMining method` -> `OUTPUT`. 

You can select one of the following combinations:

### BLAST + vote
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m vote -e 1e-03 -txl 99,97,95,90,85,80,75 -n 10 -sm 'Sample' -kp -o vote_pipe
```


### BLAST + voteSpecies
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m voteSpecies -e 1e-03 -pi 97 -n 10 -sm 'Sample' -kp -o voteSpecies_pipe
```


### BLAST + lca
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m lca -e 1e-03 -pi 97 -n 10 -sm 'Sample' -kp -o lca_pipe
```


### BLAST + besthit
```bash
$ blastMining full_pipeline -i test_data/ASV.fasta -bp "-db nt -max_target_seqs 10 -num_threads 5" -m besthit -e 1e-03 -pi 97 -n 10 -sm 'Sample' -kp -o besthit_pipe
```

---

## Command options
```bash
$ blastMining --help


usage: blastMining [-h] [-v] [--debug] {vote,voteSpecies,lca,besthit,full_pipeline} ...

blastMining v.1.0.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

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



### Method A
```bash
$ blastMining vote --help


usage: blastMining vote [-h] -i INPUT [-e EVALUE] [-txl TAXA_LEVEL] [-n TOPN] [-sm SAMPLE_NAME] [-kp] -o OUTPUT [-v]

blastMining: voting method with pident cut-off

blastMining v.1.0.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -txl TAXA_LEVEL, --taxa_level TAXA_LEVEL
                        P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
                        [default=99,97,95,90,85,80,75]
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -kp, --krona_plot     Draw krona plot
  -o OUTPUT, --output OUTPUT
                        output
  -v, --version         show program's version number and exit
```



### Method B
```bash
$ blastMining voteSpecies --help


usage: blastMining voteSpecies [-h] -i INPUT [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] [-kp] -o OUTPUT [-v]

blastMining: vote at species level for all

blastMining v.1.0.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity
                        (Ignore hits if their p. identities are below this threshold)
                        [default=97]
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -kp, --krona_plot     Draw krona plot
  -o OUTPUT, --output OUTPUT
                        output
  -v, --version         show program's version number and exit
```



### Method C
```bash
$ blastMining lca --help


usage: blastMining lca [-h] -i INPUT [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] [-kp] -o OUTPUT [-v]

blastMining: lca method

blastMining v.1.0.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity
                        (Ignore hits if their p. identities are below this threshold)
                        [default=97]
  -n TOPN, --topN TOPN  Top N hits used for LCA calculation [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -kp, --krona_plot     Draw krona plot
  -o OUTPUT, --output OUTPUT
                        output
  -v, --version         show program's version number and exit
```



### Method D
```bash
$ blastMining besthit --help


usage: blastMining besthit [-h] -i INPUT [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] [-kp] -o OUTPUT [-v]

blastMining: besthit method

blastMining v.1.0.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity
                        (Ignore hits if their p. identities are below this threshold)
                        [default=97]
  -n TOPN, --topN TOPN  Top N hits used for sorting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -kp, --krona_plot     Draw krona plot
  -o OUTPUT, --output OUTPUT
                        output
  -v, --version         show program's version number and exit
```



### Full pipeline
```bash
$ full_pipeline --help


usage: blastMining full_pipeline [-h] -i INPUT -bp BLAST_PARAM [-m MINING] [-e EVALUE] [-pi PIDENT] [-txl TAXA_LEVEL] [-n TOPN] [-sm SAMPLE_NAME] [-kp] -o OUTPUT [-v]

blastMining: Running BLAST + mining the output

blastMining v.1.0.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input FASTA
  -bp BLAST_PARAM, --blast_param BLAST_PARAM
                        BLAST parameters:
                        Note: "-outfmt" has been defined by the package, you don't need to add it
                        (-outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid")
                        [default="-db nt -num_threads 1 -max_target_seqs 10"]
  -m MINING, --mining MINING
                        blast mining method
                        Available methods={'vote','voteSpecies','lca','besthit'}
                        [default='vote']
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity
                        (Ignore hits if their p. identities are below this threshold)
                        [default=97]
                        **Not compatible** with "vote method"
  -txl TAXA_LEVEL, --taxa_level TAXA_LEVEL
                        P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
                        [default=99,97,95,90,85,80,75]
                        **Required** for "vote method"
  -n TOPN, --topN TOPN  Top N hits used for voting [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table [default="sample"]
  -kp, --krona_plot     Draw krona plot
  -o OUTPUT, --output OUTPUT
                        output
  -v, --version         show program's version number and exit
```
---
## Utility
In the case you want to convert the `OUTPUT.summary` to the `krona-input format` (OUTPUT.krona) for interactive [krona pie charts visualization](http://marbl.github.io/Krona/img/screen_phymmbl.png), 
you can use the following script to do so.

```bash
$ tab2krona.py -i OUTPUT.summary -o OUTPUT
```
The full command of the above script is as follow.
```bash
$ tab2krona.py --help


usage: tab2krona.py [-h] [-v] -i INPUT [-o OUTPUT]

convert TABLE.summary to TABLE.krona

***
This script is a part of blastMining program
***

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -v, --version         print version and exit
  -i INPUT, --input INPUT
                        input table
  -o OUTPUT, --output OUTPUT
                        output name [default = 'OUTPUT']

```
 

# Citation
**If you find this package useful**, `please cite`:
```BibTeX
@article{
 author = {Khoiri, Ahmad Nuruddin},
 title = {blastMining: mining blast output},
 year = {2022},
 DOI = {10.5281/zenodo.6823244},
 URL = { + https://github.com/NuruddinKhoiry/blastMining},
}
```
