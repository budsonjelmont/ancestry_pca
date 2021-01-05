#!/bin/bash

# Convert a .vcf file to a plink binary fileset (.bed, .bim, .fam)

module load plink/1.90

vcfdir=/sc/arion/scratch/belmoj01/splicingQTL/
vcffile=Capstone4.sel.hasPhenosOnly.idsync.vcf

mkdir -p $vcfoutdir

plink --vcf $vcfdir/$vcffile --double-id --make-bed --out $vcfdir/$vcffile
