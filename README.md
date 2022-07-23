# blastMining

### Mining BLAST OUTPUT

[![PyPI Version](https://img.shields.io/pypi/v/blastMining?style=flat-square)](https://pypi.org/project/blastMining)
[![Conda Version](https://anaconda.org/bioconda/blastmining/badges/version.svg)](https://anaconda.org/bioconda/blastmining)
[![Last updated](https://anaconda.org/bioconda/blastmining/badges/latest_release_date.svg)](https://anaconda.org/bioconda/blastmining)
[![Platform](https://anaconda.org/bioconda/blastmining/badges/platforms.svg)](https://anaconda.org/bioconda/blastmining)
[![Downloads](https://anaconda.org/bioconda/blastmining/badges/downloads.svg)](https://anaconda.org/bioconda/blastmining/files)
[![Install](https://anaconda.org/bioconda/blastmining/badges/installer/conda.svg)](https://anaconda.org/bioconda/blastmining)
[![License](https://anaconda.org/bioconda/blastmining/badges/license.svg)](https://github.com/NuruddinKhoiry/blastMining/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/504235539.svg)](https://zenodo.org/badge/latestdoi/504235539)

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

* [TaxonKit - NCBI Taxonomy Toolkit](https://bioinf.shenwei.me/taxonkit/)

* [csvtk](https://github.com/shenwei356/csvtk)

* [Python3](https://www.python.org/)

* [krona](https://github.com/marbl/Krona/wiki)

* [GNU Parallel](https://www.gnu.org/software/parallel) 


## Installation

### Option 1. Install via [conda](https://anaconda.org/bioconda/blastmining)

**NOTE:** `blastMining v.1.1.0` is not available in the conda installation yet.
 
This option will automatically install the dependecy programs. So, you don't need to install them manually. 

```bash
$ conda install -c bioconda blastmining
```
**[OPTIONAL]** In the case your `blast version` is lower than `2.12.0`, you can update it using `conda`.
```bash
$ conda install -c bioconda blast=2.12.0
```

### Option 2. Install via [PyPI](https://pypi.org/project/blastMining/)

```bash
$ pip install blastMining
```

### Option 3. Install manually

Download the latest realese of [blastMining](https://github.com/NuruddinKhoiry/blastMining/releases/download/1.1.0/blastMining-1.1.0.tar.gz) in my Github repository.

Then install it using pip

```bash
$ pip install blastMining-1.1.0.tar.gz
```

### Installation Notes

If you install `blastMining` using **option 2** or **option 3**, you need to install the dependency programs.

**You can install the dependecy programs with conda**

Make sure your conda environment is `up to date` for the sake of the dependency programs.

```bash
$ conda update -n base conda

$ conda install -c bioconda taxonkit csvtk krona blast=2.12.0

$ conda install -c conda-forge parallel
```

### Before use

**Don't forget** to install the required databases for `BLAST` and [`TaxonKit`](https://bioinf.shenwei.me/taxonkit/usage/#before-use)  


## Tutorial
Running blastn
```bash
$ blastn -query ASV.fasta \
	-db nt \
	-out BLASTn.out \
	-outfmt="6 qseqid sseqid pident length mismatch gapopen evalue bitscore staxid" \
	-max_target_seqs 10
```
**Note**: Please strict to the above `blast outfmt`

Next, `mining` your blast result with **one of the following methods**:

### * Method A. Majority vote with percent identity cut-off

The `vote` algorithm is as follow:

![vote method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/vote_method.png?raw=true)

The default percent identity cut-off is `99`, `97`, `95`, `90`, `85`, `80`, and `75` for `Species`, `Genus`, `Family`, `Order`, `Class`, `Phylum`, and `Kingdom`, respectively.

```bash
$ blastMining vote \
	-i BLASTn.out \
	-o vote_method \
	-e 0.001 \
	-txl 99,97,95,90,85,80,75 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p lca_method \
	-kp \
	-rm
```


### * Method B. Majority vote at species level

The `voteSpecies` algorithm is as follow:

![vote method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/voteSpecies_method.png?raw=true)


```bash
$ blastMining voteSpecies \
	-i BLASTn.out \
	-o voteSpecies_method \
	-e 0.001 \
	-pi 99 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p voteSpecies_method \
	-kp \
	-rm 
```


### * Method C. LCA

The `lca` algorithm in is as follow:

![lca method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/lca_method.png?raw=true) 

The `lca` algorithm used in `blastMining` is from [TaxonKit](https://bioinf.shenwei.me/taxonkit/usage/#lca).

```bash
$ blastMining lca \
	-i BLASTn.out \
	-o lca_method \
	-e 0.001 \
	-pi 95 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p lca_method \
	-kp \
	-rm 
```


### * Method D. besthit

The `besthit` algorithm is as follow:

![besthit method](https://github.com/NuruddinKhoiry/blastMining/blob/master/docs/images/besthit_method.png?raw=true) 

```bash
$ blastMining besthit \
	-i BLASTn.out \
	-o besthit_method \
	-e 0.001 \
	-pi 97 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p besthit_method \
	-kp \
	-rm
```

----

## Full_pipeline option

This option allows you to run a full pipeline started from `blastn` -> `blastn_output` ->  `blastMining method` -> `OUTPUT`. 

You can select one of the following combinations:

### BLAST + vote
```bash
$ blastMining full_pipeline \
	-i ASV.fasta \
	-o vote_pipe \
	-bp "-db nt -max_target_seqs 10 -num_threads 5" \
	-m vote \
	-e 0.001 \
	-txl 99,97,95,90,85,80,75 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p lca_method \
	-kp \
	-rm
```


### BLAST + voteSpecies
```bash
$ blastMining full_pipeline \
	-i ASV.fasta \
	-o voteSpecies_pipe \
	-bp "-db nt -max_target_seqs 10 -num_threads 5" \
	-m voteSpecies \
	-e 0.001 \
	-pi 99 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p voteSpecies_method \
	-kp \
	-rm 
```


### BLAST + lca
```bash
$ blastMining full_pipeline \
	-i ASV.fasta \
	-o lca_pipe \
	-bp "-db nt -max_target_seqs 10 -num_threads 5" \
	-m lca \
	-e 0.001 \
	-pi 99 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p lca_method \
	-kp \
	-rm 
```


### BLAST + besthit
```bash
$ blastMining full_pipeline \
	-i ASV.fasta \
	-o besthit_pipe \
	-bp "-db nt -max_target_seqs 10 -num_threads 5" \
	-m besthit \
	-e 0.001 \
	-pi 97 \
	-n 10 \
	-sm 'Sample' \
	-j 8 \
	-p besthit_method \
	-kp \
	-rm 
```

---

## Command options
```bash
$ blastMining --help


usage: blastMining [-h] [-v] {vote,voteSpecies,lca,besthit,full_pipeline} ...

blastMining v.1.1.0

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
```



### Method A
```bash
$ blastMining vote --help


usage: blastMining vote [-h] [-v] -i INPUT -o OUTDIR [-e EVALUE] [-txl TAXA_LEVEL] [-n TOPN] [-sm SAMPLE_NAME] [-j JOBS] [-p PREFIX] [-kp] [-rm]

blastMining: voting method with pident cut-off

blastMining v.1.1.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
                        [required]
  -o OUTDIR, --outdir OUTDIR
                        Output directory
                        [required]
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -txl TAXA_LEVEL, --taxa_level TAXA_LEVEL
                        P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
                        [default=99,97,95,90,85,80,75]
  -n TOPN, --topN TOPN  Top N hits used for voting
                        [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table
                        [default="sample"]
  -j JOBS, --jobs JOBS  Number of jobs to run parallelly
                        [default=1]
  -p PREFIX, --prefix PREFIX
                        Output prefix
                        [default='vote_method']
  -kp, --krona_plot     Draw krona plot
                        [default=False]
  -rm, --rm_tmpdir      Remove temporary directory (TMPDIR)
                        [default=False]
```



### Method B
```bash
$ blastMining voteSpecies --help


usage: blastMining voteSpecies [-h] [-v] -i INPUT -o OUTDIR [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] [-j JOBS] [-p PREFIX] [-kp] [-rm]

blastMining: vote at species level for all

blastMining v.1.1.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
                        [required]
  -o OUTDIR, --outdir OUTDIR
                        Output directory
                        [required]
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity
                        (Ignore hits if their p. identities are below this threshold)
                        [default=99]
  -n TOPN, --topN TOPN  Top N hits used for voting
                        [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table
                        [default="sample"]
  -j JOBS, --jobs JOBS  Number of jobs to run parallelly
                        [default=1]
  -p PREFIX, --prefix PREFIX
                        Output prefix
                        [default='voteSpecies_method']
  -kp, --krona_plot     Draw krona plot
                        [default=False]
  -rm, --rm_tmpdir      Remove temporary directory (TMPDIR)
                        [default=False]
```



### Method C
```bash
$ blastMining lca --help


usage: blastMining lca [-h] [-v] -i INPUT -o OUTDIR [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] [-j JOBS] [-p PREFIX] [-kp] [-rm]

blastMining: lca method

blastMining v.1.1.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        blast.out file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
                        [required]
  -o OUTDIR, --outdir OUTDIR
                        Output directory
                        [required]
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity
                        (Ignore hits if their p. identities are below this threshold)
                        [default=99]
  -n TOPN, --topN TOPN  Top N hits used for LCA calculation
                        [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table
                        [default="sample"]
  -j JOBS, --jobs JOBS  Number of jobs to run parallelly
                        [default=1]
  -p PREFIX, --prefix PREFIX
                        Output prefix
                        [default='lca_method']
  -kp, --krona_plot     Draw krona plot
                        [default=False]
  -rm, --rm_tmpdir      Remove temporary directory (TMPDIR)
                        [default=False]
```



### Method D
```bash
$ blastMining besthit --help


usage: blastMining besthit [-h] [-v] -i INPUT -o OUTDIR [-e EVALUE] [-pi PIDENT] [-n TOPN] [-sm SAMPLE_NAME] [-j JOBS] [-p PREFIX] [-kp] [-rm]

blastMining: besthit method

blastMining v.1.1.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        Input file. Please use this blast outfmt 6 ONLY:
                        ("qseqid","sseqid","pident","length","mismatch","gapopen","evalue","bitscore","staxid")
                        [required]
  -o OUTDIR, --outdir OUTDIR
                        Output directory
                        [required]
  -e EVALUE, --evalue EVALUE
                        Threshold of evalue
                        (Ignore hits if their evalues are above this threshold)
                        [default=1-e3]
  -pi PIDENT, --pident PIDENT
                        Threshold of p. identity
                        (Ignore hits if their p. identities are below this threshold)
                        [default=97]
  -n TOPN, --topN TOPN  Top N hits used for sorting
                        [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table
                        [default="sample"]
  -j JOBS, --jobs JOBS  Number of jobs to run parallelly
                        [default=1]
  -p PREFIX, --prefix PREFIX
                        Output prefix
                        [default='besthit_method']
  -kp, --krona_plot     Draw krona plot
                        [default=False]
  -rm, --rm_tmpdir      Remove temporary directory (TMPDIR)
                        [default=False]
```



### Full pipeline
```bash
$ full_pipeline --help


usage: blastMining full_pipeline [-h] [-v] -i INPUT -o OUTDIR -bp BLAST_PARAM [-m MINING] [-e EVALUE] [-pi PIDENT] [-txl TAXA_LEVEL] [-n TOPN] [-sm SAMPLE_NAME] [-j JOBS]
                                 [-p PREFIX] [-kp] [-rm]

blastMining: Running BLAST + mining the output

blastMining v.1.1.0

Written by: Ahmad Nuruddin Khoiri (nuruddinkhoiri34@gmail.com)

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        input FASTA
                        [required]
  -o OUTDIR, --outdir OUTDIR
                        Output directory
                        [required]
  -bp BLAST_PARAM, --blast_param BLAST_PARAM
                        BLAST parameters:
                        Note: "-outfmt" has been defined by the package, you don't need to add it
                        [default="-db nt -num_threads 1 -max_target_seqs 10"]
  -m MINING, --mining MINING
                        blastMining method
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
                        **Required** for "voteSpecies, lca, and besthit methods"
                        **Not compatible** with "vote method"
  -txl TAXA_LEVEL, --taxa_level TAXA_LEVEL
                        P.identity cut-off for Kingdom,Phylum,Class,Order,Family,Genus,Species
                        [default=99,97,95,90,85,80,75]
                        **Required** for "vote method"
                        **Not compatible** with "voteSpecies, lca, and besthit methods"
  -n TOPN, --topN TOPN  Top N hits used for voting
                        [default=10]
  -sm SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name in the print out table
                        [default="sample"]
  -j JOBS, --jobs JOBS  Number of jobs to run parallelly
                        [default=1]
  -p PREFIX, --prefix PREFIX
                        Output prefix
                        [default='blastMining']
  -kp, --krona_plot     Draw krona plot
                        [default=False]
  -rm, --rm_tmpdir      Remove temporary directory (TMPDIR)
                        [default=False]
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
 DOI = {10.5281/zenodo.6891894},
 URL = { + https://github.com/NuruddinKhoiry/blastMining},
}
```
