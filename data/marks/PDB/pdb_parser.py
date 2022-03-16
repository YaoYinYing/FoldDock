import argparse
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Extract the sequences from individual chains from a pdb file.''')

parser.add_argument('--pdbdir', nargs=1, type= str, default=sys.stdin, help = 'Path to data.')
parser.add_argument('--pdbmeta', nargs=1, type= str, default=sys.stdin, help = 'Path to pdb ids and interacting chains.')
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

def read_all_chains_seq(pdbname):
    '''Read the sequence of all chains in a pdb file
    '''
    three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
    'SEC':'U', 'PYL':'O', 'GLX':'X', 'UNK': 'X'}

    with open(pdbname) as pdbfile:
        chains_sequence = {} #Sequence
        prev_res_no=''
        for line in pdbfile:
            try:
                if not line.startswith('ATOM'):
                    continue
                atm_record = parse_atm_record(line)
                if atm_record['atm_name']!='CA' or atm_record['res_no']==prev_res_no:
                    continue
                if atm_record['chain'] in [*chains_sequence.keys()]:
                    chains_sequence[atm_record['chain']]+=three_to_one[atm_record['res_name']]
                    prev_res_no = atm_record['res_no']
                else:
                    chains_sequence[atm_record['chain']] = three_to_one[atm_record['res_name']]
                    prev_res_no = atm_record['res_no']
            except:
                pdb.set_trace()

    return chains_sequence

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


def write_fasta(merged, outdir):
    '''Write the extracted chains to fasta
    '''

    merged = merged.dropna()
    merged = merged.reset_index()
    for i in range(len(merged)):
        row = merged.loc[i]
        with open(outdir+row.Chain+'.fasta','w') as file:
            file.write('>'+row.Chain+'\n')
            file.write(row.Seq)

def write_coordinates(merged, outdir):
    '''Save the CB coordinates for later processing
    '''
    for i in range(len(merged)):
        row = merged.loc[i]
        np.save(outdir+row.Chain+'_CB.npy',row.Coordinates)



#Parse args
args = parser.parse_args()
#Data
pdbdir = args.pdbdir[0]
pdbmeta = pd.read_csv(args.pdbmeta[0])
outdir = args.outdir[0]

#Get all pdb ids and chains
ch1 = pdbmeta.PDB+'_'+pdbmeta['Chain 1']
ch2 = pdbmeta.PDB+'_'+pdbmeta['Chain 2']
all_chains = np.append(ch1.values,ch2.values)
chain_df = pd.DataFrame()
chain_df['Chain'] = all_chains

#Read all pdb files and extract the sequence for each chain
pdb_files = glob.glob(pdbdir+'*.pdb')
fetched_chains = []
fetched_seqs = []
fetched_seqlens = []
fetched_coordinates = []
for name in pdb_files:
    chains_sequence = read_all_chains_seq(name)
    chains_coordinates = read_all_chains_coords(name)
    id = (name.split('/')[-1][:4]).upper()
    for key in chains_sequence:

        fetched_chains.append(id+'_'+key)
        fetched_seqs.append(chains_sequence[key])
        fetched_seqlens.append(len(chains_sequence[key]))
        try:
            fetched_coordinates.append(np.array(chains_coordinates[key]))
        except:
            fetched_coordinates.append(np.array([]))


#Create df
fetched_df = pd.DataFrame()
fetched_df['Chain'] = fetched_chains
fetched_df['Seq'] = fetched_seqs
fetched_df['Seqlen'] = fetched_seqlens
fetched_df['Coordinates'] = fetched_coordinates
#Merge
merged = pd.merge(chain_df, fetched_df, on='Chain', how='left')

print(merged[merged.Seq.isnull()])
print(merged[merged.Coordinates.isnull()])
pdb.set_trace()
#Write fasta
#write_fasta(merged, outdir)

#Write coordinates
write_coordinates(merged, pdbdir+'CB_coordinates/')
