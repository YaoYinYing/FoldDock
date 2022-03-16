import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
from sklearn import metrics
import pdb

parser = argparse.ArgumentParser(description = '''Evaluate the DockQ scores.''')


parser.add_argument('--dockqfiles', nargs=1, type= str, default=sys.stdin, help = 'Path to dockq scores.')
parser.add_argument('--outfile', nargs=1, type= str, default=sys.stdin, help = 'Path to output file in csv')


################FUNCTIONS#################
def parse_dockq(dockqfile):
    '''Parse the DockQ scores for each complex
    '''
    names = []
    scores = []
    with open(dockqfile, 'r') as file:
        for line in file:
            if 'Native :' in line:
                names.append(line.split('/')[-1].split('.')[0])
            if 'DockQ 0.' in line:
                scores.append(float(line.split()[1]))

    #Create df
    score_df = pd.DataFrame()
    score_df['complex_id'] = names
    score_df['DockQ'] = scores
    score_df = score_df.drop_duplicates()
    return score_df

####################MAIN##############################
#Parse args
args = parser.parse_args()
#Data
dockqfiles = glob.glob(args.dockqfiles[0]+'*.txt')
outfile = args.outfile[0]

#Parse
score_df = ''
for dockqfile in dockqfiles:
    score_df_partial = parse_dockq(dockqfile)
    mode = dockqfile.split('/')[1].split('.')[0]

    #Rename
    score_df_partial = score_df_partial.rename(columns={"DockQ": "DockQ_"+mode})
    if len(score_df)>1:
        score_df = pd.merge(score_df,score_df_partial,on='complex_id',how='inner')
    else:
        score_df = score_df_partial


    print(mode,len(score_df))

#Save
score_df.to_csv(outfile,index=False)
