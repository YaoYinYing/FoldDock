import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Select matching OX identifiers for interacting chains.''')

parser.add_argument('--oxdir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with a3m matrices and OX information.')
parser.add_argument('--pdbmeta', nargs=1, type= str, default=sys.stdin, help = 'Path to pdb ids and interacting chains.')
parser.add_argument('--mode', nargs=1, type= str, default=sys.stdin, help = 'Selection mode: top/bottom/all/random100.')
parser.add_argument('--idlist', nargs=1, type= str, default=sys.stdin, help = 'A list of specific ids to use - otherwise all in pdbmeta are used.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')


def select_top_ox(ox1, ox2, msa1, msa2):
    '''Select the top ox match (first match) in each MSA and merge the sequences to a final MSA file
    in a3m format
    '''
    #Don't remove the zeros (no OX), then the query sequences (first line)
    #will be removed
    matching_ox = np.intersect1d(ox1,ox2)

    ind1 = [] #Index to select from the individual MSAs
    ind2 = []
    #Go through all matching and select the first (top) hit
    for ox in matching_ox:
        ind1.append(min(np.argwhere(ox1==ox)[:,0]))
        ind2.append(min(np.argwhere(ox2==ox)[:,0]))

    #Select from MSAs and merge
    merged = np.concatenate((msa1[ind1], msa2[ind2]),axis=1)

    return merged

def select_bottom_ox(ox1, ox2, msa1, msa2):
    '''Select the bottom ox match (last match) in each MSA and merge the sequences to a final MSA file
    in a3m format
    '''
    #Don't remove the zeros (no OX), then the query sequences (first line)
    #will be removed
    matching_ox = np.intersect1d(ox1,ox2)
    matching_ox = matching_ox[matching_ox!=0]

    ind1 = [0] #Index to select from the individual MSAs
    ind2 = [0]
    #Go through all matching and select the first (top) hit
    for ox in matching_ox:
        ind1.append(max(np.argwhere(ox1==ox)[:,0]))
        ind2.append(max(np.argwhere(ox2==ox)[:,0]))

    #Select from MSAs and merge
    merged = np.concatenate((msa1[ind1], msa2[ind2]),axis=1)

    return merged

def select_all(msa1, msa2):
    '''Select all top sequences according to the length of the shortest MSA and merge
    '''

    l = min([msa1.shape[0],msa2.shape[0]])

    #Select from MSAs and merge
    merged = np.concatenate((msa1[:l,:], msa2[:l,:]),axis=1)

    return merged

def select_random100(ox1, ox2, msa1, msa2):
    '''Select 100 random matches that are not top or bottom in each MSA and merge the sequences to a final MSA file
    in a3m format
    '''

    matching_ox = np.intersect1d(ox1,ox2)
    matching_ox = matching_ox[matching_ox!=0]

    random100 = []
    for i in range(100):
        ind1 = [0] #Index to select from the individual MSAs
        ind2 = [0]
        #Go through all matching and select the first (top) hit
        for ox in matching_ox:
            ind1.append(np.random.choice(np.argwhere(ox1==ox)[:,0],size=1, replace=False)[0])
            ind2.append(np.random.choice(np.argwhere(ox2==ox)[:,0],size=1, replace=False)[0])

        #Select from MSAs and merge
        merged = np.concatenate((msa1[ind1], msa2[ind2]),axis=1)
        random100.append(merged)

    return random100

def write_a3m(merged, outfile):
    '''Write a3m MSA'''
    backmap = { 1:'A', 2:'C', 3:'D', 4:'E', 5:'F',6:'G' ,7:'H',
               8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P',14:'Q',
               15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y', 21:'-'} #Here all unusual AAs and gaps are set to the same char (same in the GaussDCA script)

    with open(outfile,'w') as file:
        for i in range(len(merged)):
            file.write('>'+str(i)+'\n')
            file.write(''.join([backmap[ch] for ch in merged[i]])+'\n')

    return None


#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
oxdir = args.oxdir[0]
pdbmeta = pd.read_csv(args.pdbmeta[0])
mode = args.mode[0]
idlist = pd.read_csv(args.idlist[0])
outdir = args.outdir[0]

#Go through all pairs in pdbmeta

for i in range(len(pdbmeta)):
    row = pdbmeta.loc[i]
    id1 = row.PDB+'_'+row['Chain 1']
    id2 = row.PDB+'_'+row['Chain 2']
    if len(idlist)>0:
        if id1 not in idlist.id1.values and id2 not in idlist.id2.values:
            continue
    try:
        ox1 = np.load(oxdir+id1+'_species.npy')
        ox2 = np.load(oxdir+id2+'_species.npy')
        msa1 = np.load(oxdir+id1+'_a3m_matrix.npy')
        msa2 = np.load(oxdir+id2+'_a3m_matrix.npy')

    except:
        continue

    #Select
    if mode=='top':
        merged = select_top_ox(ox1, ox2, msa1, msa2)
        write_a3m(merged, outdir+id1+'_'+id2+'_top.a3m')

    if mode=='bottom':
        merged = select_bottom_ox(ox1, ox2, msa1, msa2)
        write_a3m(merged, outdir+id1+'_'+id2+'_bottom.a3m')

    if mode=='all':
        merged = select_all(msa1, msa2)
        write_a3m(merged, outdir+id1+'_'+id2+'_all.a3m')

    if mode=='random100':
        random100 = select_random100(ox1,ox2,msa1, msa2)
        for i in range(len(random100)):
            write_a3m(random100[i], outdir+id1+'_'+id2+'_'+str(i)+'.a3m')
