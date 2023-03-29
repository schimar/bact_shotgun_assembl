import os
#import numexpr
#os.environ["NUMEXPR_MAX_THREADS"]="272"
#if 'NUMEXPR_MAX_THREADS' in os.environ: os.environ.pop('NUMEXPR_MAX_THREADS')
#import numexpr
#print('NumExpr.nthreads =' + str(numexpr.nthreads))
import pandas as pd
from glob import glob
import shutil 
#import sgkit as sg
import numpy as np
#import scipy.spatial
##%matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import random
import sys

sns.set_style("whitegrid")
random.seed(44)

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info 
#samples24 = pd.read_csv(config['samplesSub'], sep='\t', index_col=False)
#
#sample_names = list(samples24['sample'])
#sample_locations = list(samples24['location'])
#samples_set = zip(sample_names, sample_locations)
#samples_dict = dict(zip(sample_names, sample_locations))
#
## load sample info for all ant samples (the 85 + 24 = 109) 
#ids = pd.read_table(config['ids'], sep='\t', index_col=False)
#ids['id_nest'] = ids['id'] + '_' + ids['nest']
#iddict = ids[['sample', 'id_nest']].set_index('sample').T.to_dict('list')



###### helper functions ######

def getFqHome(sample):
  return(list(os.path.join(samples_dict[sample],"{0}_{1}_001.fastq.gz".format(sample,pair)) for pair in ['R1','R2']))

def getTrmHome(sample):
  return(list(os.path.join('trm', "{0}_{1}.fq.gz".format(sample,pair)) for pair in ['R1','R2']))


def renameBamBai(bam, bai):
  bam = str(bam).strip('{').strip('}').strip("'")
  bamDest = 'bbmap/all/'
  bamID = bam.split('/')[1].split('.')[0]
  newBamID = iddict[bamID][0]
  shutil.copy(bam, bamDest)
  os.rename(''.join([bamDest, '.'.join([bamID, 'bam'])]), ''.join([bamDest, '.'.join([newBamID, 'bam'])]))
  # bai
  bai = str(bai).strip('{').strip('}').strip("'")
  shutil.copy(bai, bamDest)
  os.rename(''.join([bamDest, '.'.join([bamID, 'bam.bai'])]), ''.join([bamDest, '.'.join([newBamID, 'bam.bai'])]))



