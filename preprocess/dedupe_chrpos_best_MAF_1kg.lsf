#!/bin/bash

# Preprocess 1KG reference to make a reference containing only SNPs that are unique by position. When duplicates are found, the SNP with the higher MAF is kept. 

module load plink/1.90

refdir=/sc/arion/projects/pintod02c/1kg_phase3/
tmpdir=/sc/arion/scratch/belmoj01/QTL_VCF/
refname=all_phase3
uqidname=$refname.uniqueIDs
outname=$refname.dedupeByPos_bestMAF

##################################################

# If the deduplicated-by-position reference already exists, then quit without doing anything

if [ -f ${tmpdir}${outname}.bim ]; then
   echo "${tmpdir}${outname}.bim already exists--do you really want to go to all this trouble?"
   exit 1
fi

##################################################

# Make directory to do our dirty work in
mkdir -p $tmpdir

##################################################
# Copy 1kg reference to tmp folder
cp $refdir${refname}.{bed,bim,fam} $tmpdir

# Remap IDs in .bim file to unique integer IDs (row # in the .bim file)
#awk 'BEGIN{OFS=FS="\t"}; {print $1,$1":"$4":"$5":"$6":"$2"_"NR,$3,$4,$5,$6}' $tmpdir${refname}.bim > $tmpdir${refname}.newIDs.bim # Using chr:pos:allele:altallele as an ID leads to memory problems when re-mapping IDs later, so use row int instead
awk 'BEGIN{OFS=FS="\t"}; {print $1,NR,$3,$4,$5,$6}' $tmpdir${refname}.bim > $tmpdir${refname}.newIDs.bim
# Make maps of chr:pos & RSID to new integer ID 
#awk 'BEGIN{OFS=FS="\t"}; {print $1":"$4":"$5":"$6":"$2"_"NR,$2,$1":"$4}' $tmpdir${refname}.bim > $tmpdir${refname}.idmap # See note above
awk 'BEGIN{OFS=FS="\t"}; {print NR,$2,$1":"$4}' $tmpdir${refname}.bim > $tmpdir${refname}.idmap
# Swap out the old .bim file for our newly created one  & rename files
mv $tmpdir${refname}.newIDs.bim $tmpdir${uqidname}.bim
mv $tmpdir${refname}.bed $tmpdir${uqidname}.bed
mv $tmpdir${refname}.fam $tmpdir${uqidname}.fam 

# Make MAF report w/ new IDs
plink --bfile $tmpdir$uqidname \
      --allow-extra-chr \
      --freq \
      --out ${tmpdir}${uqidname}

# .frq file formatting is wonky (multiple spaces as delimiters+has leading spaces) so convert to tab-delimited first
tr -s ' ' < ${tmpdir}${uqidname}.frq | sed -e 's/^[ ]//g' | sed -e 's/[ ]/\t/g' >  ${tmpdir}${uqidname}.frq.tab

################################################
# Deduplicate chr:pos by MAF
# First, read and store the MAF for each unique SNP ID. Then read through the new ID map and store all chr:pos and SNP int IDs. When a duplicate chr:pos is found, replace the stored value with the unique ID that has the higher MAF 
awk 'BEGIN{OFS=FS="\t"};
  NR==FNR{MAFs[$2]=$5;next};
  {if($3 in chrpos){
    if(MAFs[$1] > MAFs[chrpos[$3]]){
      chrpos[$3]=$1 
    } 
   } else {
     chrpos[$3]=$1
   }
  };
  END{for(sid in chrpos){print chrpos[sid],sid}}
 ' ${tmpdir}${uqidname}.frq.tab $tmpdir${refname}.idmap > $tmpdir${refname}.idmap.deduped 

# Extract the first column of the deduped ID map since these are the plink IDs to extract from the file before we rename  
cut -d$'\t' -f2 $tmpdir${refname}.idmap.deduped > $tmpdir${refname}.idmap.deduped.extract.in

#########################################
# Filter the 1kg reference to keep just the unique IDs selected in the last step, which should now be unique by chr:pos.
# Then, remap the unique IDs to chr:pos. This leaves a reference containing 84,280,783 SNPs when using the 1KG all phase 3 data.
plink --bfile $tmpdir$uqidname \
      --allow-extra-chr \
      --extract $tmpdir${refname}.idmap.deduped.extract.in \
      --update-name $tmpdir${refname}.idmap.deduped \
      --make-bed \
      --out ${tmpdir}${outname}
