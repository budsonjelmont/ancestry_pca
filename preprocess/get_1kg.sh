#!/bin/bash

# Source: https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html

# Get 1kg phase 3 from remote source & convert to plink binary format 

refdir=/sc/arion/projects/EPIASD/splicingQTL/PCA/1kg_phase3
mkdir -p $refdir/plink_log

cd $refdir

# If DropBox links are dead, get from https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3
pgen=https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/0nz9ey756xfocjm/all_phase3.pvar.zst?dl=1
sample=https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1


# Download and decompress 1000 Genomes phase 3 data
wget $pgen
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

# Note: address doesn't work, get from here instead: https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1
wget $pvar
mv 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst

wget $sample
mv 'phase3_corrected.psam?dl=1' all_phase3.psam

# Convert 1000 Genomes phase 3 data to plink 1 binary format
ml plink2

plink2 --pfile $refdir/all_phase3 vzs\
      --max-alleles 2 \
      --make-bed \
      --out $refdir/all_phase3
mv $refdir/all_phase3.log $refdir/log
