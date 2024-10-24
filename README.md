# quantify_bam

Velocyto.py is used to quantif bam files and create three kinds of expression values for each gene: 

1. The reads that are mapping to exons only
2. The reads that map to exon/inton boundaries
3. The reads that map ambigiouse (I Need to read up on them one more time)

What this Rust program aims to do is to replicate that behaviour of the Python program, but increaes the reliability and speed of the process.

# Function

This program expects the Bam files to be sorted on chromosomal position.
It will read the Bam files feature by feature and populate a sparse matrix with <gene> <gene>_int and <gene>_amb genes measurments.

# Usage

Not finally decided at the moment. This is work in progress. Let's see if this program is written before the velocyto process has finally produced usable output.

# Missing

1. It is not compiling at the moment. Should be rather quick fix
2. The gtf library does not handle overlapping genes or strand information at the moment.
3. The end position and mapping in the gtf need to incorporate Cigar strings! ARGH!

Apart from that I hope almost everything is as it should be!