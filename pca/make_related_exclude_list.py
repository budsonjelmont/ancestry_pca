# Read in the .imiss (arg1) and .genome (arg2) files produced by plink's --missing and --genome flags respectively.
# For all rows of .imiss where IBD > 0.2 (second degree relatives), write the ID of whichever individual has the higher FMISS to SamplesToExcludeForPCA.txt

import argparse as ap
import pandas as pd
import sys
import os

parser = ap.ArgumentParser(description='Read in the .imiss (arg1) and .genome (arg2) files produced by plink\'s --missing and --genome flags respectively. For all rows of .imiss where IBD > 0.2 (second degree relatives), write the ID of whichever individual has the higher FMISS to SamplesToExcludeForPCA.txt')
parser.add_argument('imissfile', metavar='imissfile', type=str, nargs=1, help='Path to .imiss file from plink')
parser.add_argument('relfile', metavar='relfile', type=str, nargs=1, help='Path to .genome file from plink')

args = parser.parse_args()

imissfile = args.imissfile[0] #'/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/plink_reports/Capstone4.sel.idsync.maf5.mind1.geno1.hwe1e-5.highLDexcl.indep15000_150_.2.missingStats.imiss'
relfile = args.relfile[0] #'/sc/orga/projects/EPIASD/splicingQTL/intermediate_files/plink_reports/Capstone4.sel.idsync.maf5.mind1.geno1.hwe1e-5.highLDexcl.indep50_5_.5-relatedness-report.genome'

imiss = pd.read_csv(imissfile,sep='\s+')
gen = pd.read_csv(relfile,sep='\s+')

imiss.set_index('IID',inplace=True,drop=False)

rel = gen.loc[gen.PI_HAT > 0.2] 

def buildExcludeList(rows):
  fid1 = str(rows['FID1'])
  iid1 = str(rows['IID1'])
  fid2 = str(rows['FID2'])
  iid2 = str(rows['IID2'])
  fid1_fmiss = imiss.loc[iid1,'F_MISS']
  fid2_fmiss = imiss.loc[iid2,'F_MISS']
  if fid1_fmiss == fid2_fmiss:
    return ' '.join((fid1,iid1)) 
  elif fid1_fmiss > fid2_fmiss:
    return ' '.join((fid1,iid1))
  elif fid1_fmiss < fid2_fmiss:
    return ' '.join((fid2,iid2))
  else:
    print('???')

exclude = rel.apply(buildExcludeList,axis=1) 

resdir = os.path.dirname(relfile)
exclude.to_csv(resdir + '/SamplesToExcludeForPCA.txt',index=False,header=False)
