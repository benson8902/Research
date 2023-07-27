#!/bin/bash

export PATH=/media/data/scRNA_analysis/cellranger-7.1.0:$PATH

for ((i=59; i<=91; i++))
do
    SRR="SRR196879${i}"
    if [ ! -d "${SRR}" ]; then
            echo "Directory ${SRR} doesn't exist, skipping!"
            continue
    fi

    echo "Running cellranger count for ${SRR}"
    nohup cellranger count --id=run_count_pbmc_${SRR} --fastqs=./pbmc --sample=${SRR} --transcriptome=/media/data/ref/cellranger/refdata-gex-GRCh38-2020-A
done
