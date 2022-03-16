import argparse
import sys
import os
import numpy as np
import pandas as pd
from scipy.spatial import distance
from collections import Counter
from Bio import pairwise2
import pdb

parser = argparse.ArgumentParser(description = '''Get DSSP in interface.''')
parser.add_argument('--dsspdir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with DSSP results.')
parser.add_argument('--cbdir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with CB coordinates.')
parser.add_argument('--meta', nargs=1, type= str, default=sys.stdin, help = 'Path to pdb ids and interacting chains.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')

def parse_dssp(filepath):
    '''Parse DSSP
    '''

    secondary_str = {}
    residues = {}
    surface_acc = {}
    with open(filepath, 'r') as file:
        fetch_lines = False
        for line in file:
            if fetch_lines==True:
                #Fetch data
                chain = line[11]
                residue = line[13]
                ss = line[16]
                acc = line[35:38].strip()

                #Save
                if chain in [*secondary_str.keys()]:
                    secondary_str[chain].append(ss)
                    residues[chain].append(residue)
                    surface_acc[chain].append(acc)
                else: #Assign new chain
                    secondary_str[chain]=[ss]
                    residues[chain]=[residue]
                    surface_acc[chain]=[acc]

            if '#' in line:
                fetch_lines = True


    return secondary_str, residues, surface_acc

def realign(org_seq, dssp_seq):
    '''Realign the sequence extracted from DSSP to match that of the original sequence
    '''

    aln = pairwise2.align.globalxx(org_seq, dssp_seq)

    #Return the positions that are not gaps in the original aln
    return np.argwhere(np.array([*aln[0][0]])!='-')[:,0]



def get_interface_contacts(cb1, cb2, t):
    '''Get the CBs that are within threshold t Ångström
    '''

    contact_res1 = []
    contact_res2 = []
    #The contacts will be in lower triangular fashion
    for i in range(len(cb1)):
        for j in range(len(cb2)):
            dist = distance.euclidean(cb1[i], cb2[j])
            if dist < t:
                contact_res1.append(i)
                contact_res2.append(j)

    return np.unique(contact_res1), np.unique(contact_res2)



#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
dsspdir = args.dsspdir[0]
cbdir = args.cbdir[0]
meta = pd.read_csv(args.meta[0])
outdir = args.outdir[0]
t=8 #Threshold in Ångström

#Go through all pairs in meta
results = {
'id1': [],
'id2': [],
'G': [], #3-turn helix (310 helix). Min length 3 residues.
'H': [], #4-turn helix (α helix). Minimum length 4 residues.
'I': [], #5-turn helix (π helix). Minimum length 5 residues.
'T': [], #hydrogen bonded turn (3, 4 or 5 turn)
'E': [], #extended strand in parallel and/or anti-parallel β-sheet conformation. Min length 2 residues.
'B': [], #residue in isolated β-bridge (single pair β-sheet hydrogen bond formation)
'S': [], #bend (the only non-hydrogen-bond based assignment).
'C': [], #coil (residues which are not in any of the above conformations).
' ': [],
'num_if_res':[]
}

for i in range(len(meta)):
    print(i,'out of',len(meta))
    row = meta.loc[i]

    #Get pdb id
    pdb_id = row.id1.split('_')[0]
    #Get chains
    ch1 = row.id1.split('_')[1]
    ch2 = row.id2.split('_')[1]
    #Read dssp
    secondary_str, residues, surface_acc = parse_dssp(dsspdir+pdb_id+'.pdb.dssp')
    #Get dssp for chains
    ss1, res1, acc1 = np.array(secondary_str[ch1]), np.array(residues[ch1]), np.array(surface_acc[ch1])
    ss2, res2, acc2 = np.array(secondary_str[ch2]), np.array(residues[ch2]), np.array(surface_acc[ch2])
    #Get CBs
    cb1 = np.load(cbdir+row.id1+'_CB.npy')
    cb2 = np.load(cbdir+row.id2+'_CB.npy')
    #Check lens
    if len(ss1)<len(cb1) or len(ss2)<len(cb2):
        print('Missing ss assignment', row.id1, row.id2)
        continue
    try:
        if len(ss1)>len(cb1):
            #Realign
            keep_pos = realign(row.seq1,''.join(res1))
            ss1, res1, acc1 = ss1[keep_pos], res1[keep_pos], acc1[keep_pos]
            print(len(ss1),len(cb1))
        if len(ss2)>len(cb2):
            #Realign
            keep_pos = realign(row.seq2,''.join(res2))
            ss2, res2, acc2 = ss2[keep_pos], res2[keep_pos], acc2[keep_pos]
            print(len(ss2),len(cb2))
    except:
        print('Missing ss assignment', row.id1, row.id2) 
        continue


    #Get interface contacts
    contact_res1, contact_res2 = get_interface_contacts(cb1, cb2, t)
    if_ss1 = ss1[contact_res1]
    if_ss2 = ss2[contact_res2]
    #Count ss in if
    if_counts = Counter(np.concatenate([if_ss1,if_ss2]))
    total_res = len(if_ss1)+len(if_ss2)
    results['num_if_res'].append(total_res)
    #Save
    results['id1'].append(row.id1)
    results['id2'].append(row.id2)
    for key in ['G','H','I','T','E','B','S','C',' ']:
        if key in [*if_counts.keys()]:
            results[key].append(if_counts[key]/total_res)
        else:
            results[key].append(0)

#Create df
results_df = pd.DataFrame.from_dict(results)
results_df.to_csv(outdir+'dssp_df.csv')
pdb.set_trace()
