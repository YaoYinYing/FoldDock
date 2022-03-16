import argparse
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Rewrite the AF2 output to contain two different chains in the PDB files.''')

parser.add_argument('--pdbdir', nargs=1, type= str, default=sys.stdin, help = 'Path to data.')
parser.add_argument('--meta', nargs=1, type= str, default=sys.stdin, help = 'Path to pdb ids and interacting chains.')
parser.add_argument('--mode', nargs=1, type= str, default=sys.stdin, help = 'Mode: bench4 or marks.')
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
    '''Get all atom coordinates for all chains
    '''

    with open(pdbname) as pdbfile:
        pdb_chains = {} #Coordinates
        prev_res_no=''
        prev_atm = ''
        for line in pdbfile:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            if record['atm_name']==prev_atm and record['res_no']==prev_res_no:
                continue
            if record['chain'] in [*pdb_chains.keys()]:
                pdb_chains[record['chain']].append(line)
                prev_res_no= record['res_no']
                prev_atm = record['atm_name']
            else:
                pdb_chains[record['chain']] = [line]
                prev_res_no= record['res_no']
                prev_atm = record['atm_name']


    return pdb_chains

def write_pdb(chains, chain_break, outname):
    '''Save the CB coordinates for later processing
    '''


    with open(outname,'w') as file:
        #Write the first pdb file
        atomc=1
        aac=1
        chain = 'A'
        prev_res_no = 1
        for line in chains:
            record = parse_atm_record(line)
            #Update chain and reset aac
            if aac>chain_break and chain!='B':
                chain='B'
                aac=1

            #Check res number
            if record['res_no']>prev_res_no:
                aac+=1

            atom_blank = ' '*(5-len(str(atomc)))
            res_blank = ' '*(4-len(str(aac)))
            outline = line[:6]+atom_blank+str(atomc)+line[11:21]+chain+res_blank+str(aac)+line[26:]
            file.write(outline)

            prev_res_no = record['res_no']

            #Next atom
            atomc+=1




################MAIN###############
#Parse args
args = parser.parse_args()
#Data
pdbdir = args.pdbdir[0]
meta = pd.read_csv(args.meta[0])
mode = args.mode[0]
outdir = args.outdir[0]

#Read an rewrite all pdb files
if mode=='bench4':
    for i in range(len(meta)):
        row = meta.loc[i]
        #Check if file is present
        files = glob.glob(pdbdir+row.PDB+'_u1-'+row.PDB+'_u2/*.pdb')
        if len(files)>0:
            for pdbname in files:
                chains = read_all_chains_coords(pdbname)
                #Rewrite the files
                write_pdb(chains['A'], row['Sequence Length 1'], outdir+row.PDB+'_u1-'+row.PDB+'_u2.pdb')

if mode=='marks':
    for i in range(len(meta)):
        row = meta.loc[i]
        #Check if file is present
        files = glob.glob(pdbdir+row.id1+'-'+row.id2+'/*.pdb')
        if len(files)>0:
            for pdbname in files:
                chains = read_all_chains_coords(pdbname)
                #Rewrite the files
                write_pdb(chains['A'], row.l1, outdir+row.id1+'-'+row.id2+'.pdb')

if mode=='new_dimers':
    for i in range(len(meta)):
        row = meta.loc[i]
        #Check if file is present
        name = row.PDB+'_'+row['Chain 1']+'-'+row.PDB+'_'+row['Chain 2']
        files = glob.glob(pdbdir+name+'/*.pdb')
        if len(files)>0:
            for pdbname in files:
                chains = read_all_chains_coords(pdbname)
                #Rewrite the files
                write_pdb(chains['A'], row.l1, outdir+name+'.pdb')
