import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Fetch the plDDT and calculate several metrics across each chain, both chains and the interface.''')


parser.add_argument('--plDDTdir', nargs=1, type= str, default=sys.stdin, help = 'plDDTs to merge.')
parser.add_argument('--outname', nargs=1, type= str, default=sys.stdin, help = 'Path to outname')

#################MAIN####################

#Parse args
args = parser.parse_args()
outname = args.outname[0]
#Get plDDT
plDDTs = glob.glob(args.plDDTdir[0]+'*.csv')

plddt_df = ''
for i in range(len(plDDTs)):
    df = pd.read_csv(plDDTs[i])
    df = df.drop(columns=['Unnamed: 0']) #Drop index
    #Rename columns
    for col in df.columns[2:-2]:
        df = df.rename(columns={col: col+'_'+plDDTs[i].split('_')[-1][0]})

    df = df[df.columns[:-2]]
    if len(plddt_df)>0:
        plddt_df = pd.merge(plddt_df,df,on=['id1','id2'],how='inner')
    else:
        plddt_df = df

#Save
pdb.set_trace()
plddt_df.to_csv(outname)
