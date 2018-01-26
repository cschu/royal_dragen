#!/bin/bash -e
source vcftools-0.1.13

vcfin=$(sed -n "$SLURM_ARRAY_TASK_ID"p $1)
vcfout=$(dirname $vcfin)/$(basename $vcfin .gz).nod.gz
gzip -dc $vcfin | grep -v '##DRAGEN' | bgzip -c > $vcfout && tabix -p vcf $vcfout;

