#!/bin/bash

# Merge plink formatted target dataset file w/ 1kG reference (also in plink binary fileset format)
# Remove linked sited first, then QC both datasets. Align SNPs common to both datasets, & check for chr/pos mismatches & allele flips. Also remove related individuals (pi_hat > 0.2) in the reference (you should've done this already on the target dataset though). 
# Based on vignette here: https://meyer-dfflab-cshl.github.io/plinkQC/articles/AncestryCheck.html
ml plink
ml python/3.7.3 # Used for py script to compare related individuals & generate list of relatives to drop

datdir=/sc/arion/projects/EPIASD/splicingQTL/output/geno_wasp/geno_PCA
name=Capstone4.sel.idsync.2allele.maf01.mind05.geno05.hwe1e-6.deduped.COPY
refdir=/sc/arion/scratch/belmoj01/QTL_VCF
refname=all_phase3.dedupeByPos_bestMAF
reference=1kg_phase3_bestMAF
highld=/sc/hydra/projects/pintod02c/reference-databases/high_LD_regions/high_ld_and_autsomal_regions_hg19.txt
popfile=/sc/hydra/projects/pintod02c/1kg_phase3/1kg_phase3_samplesuperpopinferreddata-FID0.txt

pyscrpath=/sc/arion/projects/EPIASD/ancestry_pca/pca

# QC parameters
ldwindow=50
ldstep=5
ldr2=0.2
maf=0.05
mind=0.05
geno=0.05

####################################################################################

mkdir -p $datdir/plink_log

# Prune study data by pruning sites in LD & also removing pre-computed high-LD areas
plink --bfile $datdir/$name \
      --exclude range $highld \
      --indep-pairwise $ldwindow $ldstep $ldr2 \
      --out $datdir/$name
mv  $datdir/$name.prune.log $datdir/plink_log/$name.prune

plink --bfile $datdir/$name \
      --autosome \
      --biallelic-only --maf $maf --mind $mind \
      --extract $datdir/$name.prune.in \
      --make-bed \
      --out $datdir/$name.pruned
mv  $datdir/$name.pruned.log $datdir/plink_log/$name.pruned

# Filter reference data for the same SNP set as in study
plink --bfile $refdir/$refname \
      --allow-extra-chr \
      --biallelic-only --maf $maf --mind $mind \
      --extract $datdir/$name.prune.in \
      --make-bed \
      --out $datdir/$refname.pruned
mv  $datdir/$refname.pruned.log $datdir/plink_log/$refname.pruned

# Check and correct chromosome mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $datdir/$name.pruned.bim $datdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $datdir/$refname.toUpdateChr

plink --bfile $datdir/$refname.pruned \
      --update-chr $datdir/$refname.toUpdateChr 1 2 \
      --make-bed \
      --out $datdir/$refname.updateChr
mv $datdir/$refname.updateChr.log $datdir/plink_log/$refname.updateChr.log

# Find variants with position mismatches
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $datdir/$name.pruned.bim $datdir/$refname.pruned.bim > \
    $datdir/${refname}.toUpdatePos

# Find allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $datdir/$name.pruned.bim $datdir/$refname.pruned.bim > \
    $datdir/$refname.toFlip

# Update positions and flip alleles
plink --bfile $datdir/$refname.updateChr \
      --update-map $datdir/$refname.toUpdatePos 1 2 \
      --flip $datdir/$refname.toFlip \
      --make-bed \
      --out $datdir/$refname.flipped
mv $datdir/$refname.flipped.log $datdir/plink_log/$refname.flipped.log

# Find related individuals in the reference and remove the one with lower genotyping rate
plink --bfile $datdir/$refname.flipped --missing --out  $datdir/$refname.flipped.missingStats
plink --bfile $datdir/$refname.flipped --genome --out $datdir/$refname.flipped.relatedness
# Make list of people to drop (they'll be dropped in the next step when we remove mismatches)
python3 $pyscrpath/make_related_exclude_list.py $datdir/$refname.flipped.missingStats.imiss $datdir/$refname.flipped.relatedness.genome

# Remove mismatches
# Any alleles that do not match after allele flipping, are identified and removed from the reference dataset.
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $datdir/$name.pruned.bim $datdir/$refname.flipped.bim > \
    $datdir/$refname.mismatch

plink --bfile $datdir/$refname.flipped \
      --remove $datdir/SamplesToExcludeForPCA.txt \
      --exclude $datdir/$refname.mismatch \
      --make-bed \
      --out $datdir/$refname.clean
mv $datdir/$refname.clean.log $datdir/plink_log/

# Merge study genotypes and reference data
# The matching study and reference dataset can now be merged into a combined dataset with plink –bmerge. If all steps outlined above were conducted successfully, no mismatch errors should occur.
# After merger filter genotyping rate again 

plink --bfile $datdir/$name.pruned  \
      --bmerge $datdir/$refname.clean.bed $datdir/$refname.clean.bim \
         $datdir/$refname.clean.fam  \
      --geno $geno \
      --make-bed \
      --out $datdir/$name.merge.$refname
mv $datdir/$name.merge.$refname.log $datdir/plink_log/

# PCA on the merged data
# We can now run principal component analysis on the combined dataset using plink –pca which returns a .eigenvec file with the family and individual ID in columns 1 and 2, followed by the first 20 principal components.

plink --bfile $datdir/$name.merge.$refname \
      --pca \
      --genome \
      --out $datdir/$name.$reference #--within $popfile --pca-cluster-names AFR SAS EAS EUR AMR # See section Dimension reduction here: https://www.cog-genomics.org/plink/1.9/strat
mv $datdir/$name.$reference.log $datdir/plink_log

#plink --bfile $datdir/$name.merge.$refname \
#      --read-genome $datdir/$name.$reference.genome \
#      --cluster \
#      --ppc 1e-3 \
#      --mds-plot 2 \
#      --out $datdir/$name.$reference
