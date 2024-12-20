[![Rust](https://github.com/stela2502/quantify_bam/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/stela2502/quantify_bam/actions/workflows/rust.yml)

# quantify_bam

Velocyto.py is used to quantif bam files and create three kinds of expression values for each gene: 

1. The reads that are mapping to exons only
2. The reads that map to exon/inton boundaries
3. The reads that map ambigiouse (I Need to read up on them one more time)

What this Rust program aims to do is to replicate that behaviour of the Python program, but increases the reliability and speed of the process.

# Function

This program expects the Bam files to be sorted on chromosomal position.
It will then copare each bam match to the provided gtf information. Reads matching to the gene/exon entries will be repoirted as \<gene name>,
Reads matching to exon+intron or only intron will be reported as \<gene name>\_unspliced.
If the only matching region in the gtf is an exon, but the read is actally expressed in antisense the read will be reported as\<\_antiense> 
and is that antisense reads overlapps with the exon intron boundary it is reported as antisense_unspliced. Might not be the most clever way to treat this.

I have developed this tool espicially because Velocyto did loose the Igkc expression in our dataset. We found out that the most likely cause to this that Vecloyto classed these reads as ambigiouse as ther is a GM overlapping with the Igkc transcript. Igkc lies in an intron of Gm30211 and it's expression if therefore likely classed as ambigiouse. Here ambigiouse mapping is unhandled at the moment as this tool will always give a match to an exon before an intron and I know of no exon that is shared by two genes.

The UMI checked single cell epression matrix is reported as 10x formated triplet of files: matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz.
The different read classes have to be split apart before analyzing the data. An example of how to do that is in this [Postprocessing Notebook](notebooks/postprocessing_example.ipynb).


# Installation

You need the Rust compiler: https://www.rust-lang.org/tools/install


Then you can clone this repo and complie the code (example for a Linux system).
But it also compiles on Windows. I just never use that for actual work.

Afterwards you can use the cargo tool to compile and install the program:

```
cargo install  --git https://github.com/stela2502/quantify_bam
```

# Usage

```
quantify_gtf -h
quantify_bam 0.4.0
Stefan L. <stefan.lang@med.lu.se>

USAGE:
    quantify_gtf [OPTIONS] --bam <BAM> --gtf <GTF> --outpath <OUTPATH> --min-umi <MIN_UMI>

OPTIONS:
    -b, --bam <BAM>              the bam file to quantify
    -g, --gtf <GTF>              the gtf file fitting to the Bam file
    -h, --help                   Print help information
    -m, --min-umi <MIN_UMI>      used processor cores (default all)
    -n, --num-proc <NUM_PROC>    the minimum (UMI) reads per cell (sample + genes + antibody
                                 combined)
    -o, --outpath <OUTPATH>      the outpath
    -V, --version                Print version information
```

Not finally decided at the moment. This is work in progress. Let's see if this program is written before the velocyto process has finally produced usable output.


# Speed

Velocypt.py did take an extended amount of time to process my bam file (>3.10^9 reads). Far longer than one day on our analysis system (CPU 2 x AMD 7413 (2.65 Ghz, 2 x 24-core; Memory  256 GB). quantify_gtf did analyze this data in less than 5 hours on that same system.

