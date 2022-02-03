## README

In human genetics, an important task is validating mutations that have been found using exome or whole genome sequencing with a PCR spanning the abnormal region. This task is not difficult, but surprisingly laborious when done manually.

This workflow automates the design of PCR and qPCR primers. As input, it takes either:

- an HGVS-formatted mutation, e.g. `NM_000546.6:c.215C>G` (for PCR)
- a gene name and exon number, e.g. `NM_007294.4:14` (for qPCR)

The workflow will:

- validate [HGVS](https://varnomen.hgvs.org/bg-material/simple/) syntax
- for mutations in coding regions try to span the exon
- constrain primers to desired parameters
- construct a second backup pair of primers that does not overlap the first pair
- avoid placing primers across known common variants
- validate all primers to not produce unspecific amplicons


## Setup

```bash
conda install -y -n primer -c bioconda nextflow && conda activate primer
git clone https://github.com/phiweger/primer && cd primer
mkdir data
# ^ add data in there, see below
```

Then set the absolute (!) path to the data in `nextflow.config`. If you use other reference genomes/ annotations/ variant databases, you also have to adjust the filenames in `config.json`.



## Run

```bash
cat input.csv
# name,method,variant
# ABCA4_v1,PCR,NM_000350.3:c.4234C>T
# ABCA4_v2,PCR,NM_000350.3:c.4773+3A>G

nextflow run workflow/main.nf --input input.csv --results designs
```


## Data provenance

You need to get three things, which need to correspond to one another:

1. reference genome
2. SNPs
3. annotation

Below, we'll use the files corresponding to the human genome `hg19` (v13, [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/), last access 2021-02-03).


```bash
# 1. Reference genome and index for random access:
# - GRCh37_latest_genomic.fna.gz
# - GRCh37_latest_genomic.fna.gz.fai

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
samtools faidx GRCh37_latest_genomic.fna.gz
```


```bash
# 2. Variants and index:
# - GRCh37_latest_dbSNP_all.vcf.gz
# - GRCh37_latest_dbSNP_all.vcf.gz.tbi

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz
gunip GRCh37_latest_dbSNP_all.vcf.gz
bgzip -c GRCh37_latest_dbSNP_all.vcf > GRCh37_latest_dbSNP_all.vcf.gz
tabix -p vcf GRCh37_latest_dbSNP_all.vcf.gz
```


```bash
# 3. Annotation database:
# - hg19-p13_annotation.db
```

Get it from [OSF](https://osf.io) (project ID [7csav](https://osf.io/7csav/)); here is how we created it:


```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
gunzip GRCh37_latest_genomic.gff.gz
# switch to Python
```

```python
import gffutils
# http://daler.github.io/gffutils/
# https://daler.github.io/gffutils/database-import.html

fp = 'GRCh37_latest_genomic.gff'
db = gffutils.create_db(
    fp, 
    dbfn='hg19-p13_annotation.db',
    force=True,
    keep_order=True,
    merge_strategy='merge',
    sort_attribute_values=True)
# Takes a couple of hours on a regular laptop, then use like
db = gffutils.FeatureDB('hg19-p13_annotation.db', keep_order=True)
db['exon-NR_024540.1-3']
```


## Design choices

For future reference, below are documented some peculiarities:

- `gffutils.FeatureDB` will not accept symlinks; this means we cannot pass the database in a `nextflow` channel. 
- There are two files with parameters, one for `nextflow` (`nextflow.config`) and one to guide the primer design (`config.json`). There are more elegant solutions, but meh.
