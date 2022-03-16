import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
from collections import Counter
import pdb

parser = argparse.ArgumentParser(description = '''Calculate Neff.''')

parser.add_argument('--a3mdir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with merged a3m files.')
parser.add_argument('--pdbmeta', nargs=1, type= str, default=sys.stdin, help = 'Path to pdb ids and interacting chains.')
parser.add_argument('--mode', nargs=1, type= str, default=sys.stdin, help = 'Selection mode: top/bottom/all/random100.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')

def read_a3m(infile):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
             'G': 6,'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,'N': 12,
             'O': 21, 'P': 13,'Q': 14, 'R': 15, 'S': 16, 'T': 17,
             'V': 18, 'W': 19, 'Y': 20,'U': 21, 'Z': 21, 'X': 21, 'J': 21}

    parsed = []#Save extracted msa

    with open(infile, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('>'): #OX=OrganismIdentifier
                continue
            line = line.rstrip()

            parsed.append([mapping.get(ch, 22) for ch in line if not ch.islower()])

    return np.array(parsed, dtype=np.int8, order='F')

def calc_neff(msa, t):
    '''Calculate Neff using cutoff t
    Neff is essentially equal to the
    number of non-redundant sequences (sequence identity<0.8) in the MSA normalized
    by the query length
    https://storage.googleapis.com/plos-corpus-prod/10.1371/journal.pcbi.1007411/2/pcbi.1007411.s001.pdf?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=wombat-sa%40plos-prod.iam.gserviceaccount.com%2F20210816%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20210816T070710Z&X-Goog-Expires=86400&X-Goog-SignedHeaders=host&X-Goog-Signature=a68350f11d44be91f3de5348ffd024e739a29b6531837c8a55778d8bb28e19340880187cc6a9479bd7cfd1742936bcd7288d768e07e9cc751181f736fde530a392daee6889fcca2d9d5e7acbdb78c47beb14c8d9b8f4a0befa72435d56be51ce149277552216a4d9f0eb02795ad888e74be8ccb9426ccbd0f18fd1e1aa4c59115c240467389694fe2f4ebecb9b1bdca63e3c3c9fe2f5877bc71e063f9af24c8deb5d2ffe1212463020f06f245cf851b954be8e39a003b23bafa56babc656a16c44beeeddc3cbb05a289a4c92eca13ba95fb2d4d64d5f2bacf68be73f7ede5bda044d30ae2b4c6999dc7faf6a6821ed0e977e80aab0ec691190ed14d8c51611cd
    '''

    l = msa.shape[1]
    lt = int(l*t) #The number of positions that may differ in the cluster
    norm = 1/np.sqrt(l)

    #Cluster
    n_clusters = 0
    clustered_seqs = 0
    remaining_msa = msa
    #Go through all seqs and cluster
    while remaining_msa.shape[0]>1:
        msa_diff = remaining_msa-remaining_msa[0]
        sim_counts = Counter(np.argwhere(msa_diff==0)[:,0]) #Count the similar positions
        vals = np.array([*sim_counts.values()])
        keys = np.array([*sim_counts.keys()])
        #Get cluster
        cluster = keys[np.argwhere(vals>=lt)] #The cluster contains all sequences with over lt similarity
        n_clusters += 1
        if cluster.shape[0]==0:
            clustered_seqs +=1
            remaining_msa = remaining_msa[1:,:]
        else:
            clustered_seqs += cluster.shape[0]
            sel_inds = np.zeros(msa_diff.shape[0])
            sel_inds[cluster[:,0]]=1
            #Update remaining seqs
            remaining_msa = remaining_msa[np.argwhere(sel_inds==0)[:,0]]
        print(n_clusters, clustered_seqs)

    #Calc Neff - add the last seq
    if len(remaining_msa)>0:
        n_clusters+=1
    neff = norm*n_clusters

    return neff


#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
a3mdir = args.a3mdir[0]
pdbmeta = pd.read_csv(args.pdbmeta[0])
mode = args.mode[0]
outdir = args.outdir[0]

#Similarity cutoff
t=0.8

#Go through all pairs in pdbmeta
ids1 = []
ids2 = []
neffs = []
for i in range(len(pdbmeta)):
    row = pdbmeta.loc[i]
    id1 = row.PDB+'_'+row['Chain 1']
    id2 = row.PDB+'_'+row['Chain 2']

    try:
        msa = read_a3m(a3mdir+id1+'_'+id2+'_'+mode+'.a3m')

    except:
        continue


    neff = calc_neff(msa, t)
    #Save
    ids1.append(id1)
    ids2.append(id2)
    neffs.append(neff)
    print(neff)

#Create df
results_df = pd.DataFrame()
results_df['id1']=ids1
results_df['id2']=ids2
results_df['Neff']=neffs
results_df.to_csv(outdir+'neff_'+mode+'.csv')
