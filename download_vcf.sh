#!/bin/bash

MY_ID="your ID"
TARGET_DIR="/hpctmp/${MY_ID}/5004/data/vcf"
mkdir -p $TARGET_DIR
cd $TARGET_DIR

for i in {1..22}; do
    echo "Downloading Chromosome ${i}..."
    wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
done