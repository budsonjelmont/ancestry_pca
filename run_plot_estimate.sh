#!/bin/bash

##### User-supplied params ######
datdir=/sc/arion/scratch/belmoj01/pca_test/
datfile=Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY
refdir=/sc/arion/scratch/belmoj01/QTL_VCF
reffile=all_phase3.dedupeByPos_bestMAF
highld=/sc/hydra/projects/pintod02c/reference-databases/high_LD_regions/high_ld_and_autsomal_regions_hg19.txt
kgpopfile=/sc/hydra/projects/pintod02c/1kg_phase3/1kg_phase3_samplesuperpopinferreddata-FID0.txt

pcafile=$datfile.$reffile

centerscalebool=true

pyscrpath=/sc/arion/projects/EPIASD/ancestry_pca/pca

if [ $centerscalebool == 'true' ]
then
  centerscale='--centerscale'
else
  centerscale=''
fi

#################################

# Do PCA
pca/merge_1kg_geno_pca.sh $datdir $datfile $refdir $reffile $highld $kgpopfile $pyscrpath $pcafile

# Plot results
Rscript viz/viz_1kg_geno_pca.R $datdir $pcafile $kgpopfile --conflevel 0.95

# Estimate ancestry
Rscript est_ancestry/infer_from_1kg.R $datdir $pcafile $kgpopfile $centerscale
