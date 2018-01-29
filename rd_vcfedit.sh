#!/bin/bash -e
source vcftools-0.1.13

# uncomment next (and comment next next) line to switch between slurm array and single job
# vcfin=$(sed -n "$SLURM_ARRAY_TASK_ID"p $1)
vcfin=$1

vcfout=$(dirname $vcfin)/$(basename $vcfin .gz).nod.gz
gzip -dc $vcfin | grep -v '##DRAGEN' | bgzip -c > $vcfout && tabix -p vcf $vcfout;

