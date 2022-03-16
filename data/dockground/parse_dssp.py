import argparse
import sys
import os
import numpy as np
import pandas as pd
from scipy.spatial import distance
from collections import Counter
from Bio import pairwise2
from collections import defaultdict
import pdb

parser = argparse.ArgumentParser(description = '''Get DSSP in interface.''')
parser.add_argument('--dsspdir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with DSSP results.')
parser.add_argument('--pdbdir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with PDB files.')
parser.add_argument('--meta', nargs=1, type= str, default=sys.stdin, help = 'Path to pdb ids and interacting chains.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')


def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record

def read_all_chains_coords(pdbname):
    '''Get all atom coordinates for the CBs (CA for GLY)
    '''

    with open(pdbname) as pdbfile:
        chains_coordinates = {} #Coordinates
        prev_res_no=''
        for line in pdbfile:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            if record['res_no'] == prev_res_no:
                continue

            if record['atm_name'] == 'CB':
                pos = np.array([record['x'], record['y'], record['z']])
            elif record['atm_name'] == 'CA' and record['res_name'] == 'GLY':
                pos = np.array([record['x'], record['y'], record['z']])
            else:
                continue

            if record['chain'] in [*chains_coordinates.keys()]:
                chains_coordinates[record['chain']].append(np.array([record['x'], record['y'], record['z']]))
                prev_res_no= record['res_no']
            else:
                chains_coordinates[record['chain']] = [np.array([record['x'], record['y'], record['z']])]
                prev_res_no= record['res_no']


    return chains_coordinates


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

    return contact_res1, contact_res2



#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
dsspdir = args.dsspdir[0]
pdbdir = args.pdbdir[0]
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
'num_if_res':[],
'num_if_contacts_total':[]
}

for i in range(len(meta)):
    print(i,'out of',len(meta))
    row = meta.loc[i]
    #Read dssp
    secondary_str, residues, surface_acc = parse_dssp(dsspdir+row.PDB+'_bc.pdb.dssp')
    #Get dssp for chains
    ss1, res1, acc1 = np.array(secondary_str['A']), np.array(residues['A']), np.array(surface_acc['A'])
    ss2, res2, acc2 = np.array(secondary_str['B']), np.array(residues['B']), np.array(surface_acc['B'])
    #Get CBs
    cb_coords = read_all_chains_coords(pdbdir+row.PDB+'_bc.pdb')
    cb1, cb2 = np.array(cb_coords['A']), np.array(cb_coords['B'])

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
    if_ss1 = ss1[np.unique(contact_res1)]
    if_ss2 = ss2[np.unique(contact_res2)]
    #Count ss in if
    if_counts = Counter(np.concatenate([if_ss1,if_ss2]))
    total_res = len(if_ss1)+len(if_ss2)
    results['num_if_res'].append(total_res)
    results['num_if_contacts_total'].append(len(contact_res1)+len(contact_res2))
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
