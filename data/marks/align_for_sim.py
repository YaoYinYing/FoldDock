import argparse
import sys
import os
import numpy as np
import pandas as pd
from Bio import pairwise2
import pdb

parser = argparse.ArgumentParser(description = '''Align sequence 1 to sequence 2 to check the similarity between them.''')
parser.add_argument('--meta', nargs=1, type= str, default=sys.stdin, help = 'Path to ids and sequences.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')

#################MAIN####################

#Parse args
args = parser.parse_args()
meta = pd.read_csv(args.meta[0])
outdir = args.outdir[0]

#Go through all pairs in meta
scores = []
for i in range(len(meta)):
    print(i,'out of',len(meta))
    row = meta.loc[i]

    #Align
    scores.append(pairwise2.align.globalxx(row.seq1,row.seq2,score_only=True))
#Update df
meta['aln_score'] = scores
meta.to_csv(args.meta[0])
