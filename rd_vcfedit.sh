#!/bin/bash -e
source vcftools-0.1.13

# uncomment next (and comment next next) line to switch between slurm array and single job
# vcfin=$(sed -n "$SLURM_ARRAY_TASK_ID"p $1)
vcfin=$1

vcfout=$(dirname $vcfin)/$(basename $vcfin .gz).nod.gz
gzip -dc $vcfin | awk '/^#/ { if (NR==1 || substr($0,2,1) != "#") print $0; next;} {print $0;}' | bgzip -c > $vcfout && tabix -p vcf $vcfout;

