#!/bin/bash
#SBATCH --ntasks-per-node 48
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH -A csens2024-3-2
#SBATCH -J veloProB
#SBATCH -o veloProB.%j.out
#SBATCH -e veloProB.%j.err


samples=("ProB" "CLP" "FrA")
this=${samples[0]}

opath="${this}/velocity_benchy/"
ifile="./${this}/outs/possorted_genome_bam.bam"
gtf="/scale/gr01/shared/common/genome/10Xindexes/cellranger/9.0/mus_muculus/GRCm39/data/reference_sources/gencode.vM36.annotation.gtf.gz"

cmd="velocyto run \
        --samtools-threads 48 \
        --samtools-memory 5000 \
        -o ${opath} \
        ${ifile} \
        ${gtf} \
"

if [ ! -d $opath ]; then
        mkdir -p $opath
fi

echo ${cmd} > ${path}/start.time
date >> ${path}/start.time

echo $cmd
apptainer exec -B /scale /scale/gr01/shared/common/software/ScanpyEcoSys/1.0/ScanpyEcoSys_v1.0.sif $cmd

echo ${cmd} > ${path}/finish.time
date >> ${path}/finish.time

exit 0
