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


def read_all_chains(pdbname):
    '''Get all atom coordinates
    '''

    with open(pdbname) as pdbfile:
        chains = {} #Coordinates
        prev_res_no = ''
        prev_atom = ''
        for line in pdbfile:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            if record['res_no'] == prev_res_no and record['atm_name']==prev_atom:
                continue

            if record['atm_name']=='CB':
                #Move the previous atom (O)
                #AF2 writes N, CA, C, CB, O although PDB is N, CA, C, O, CB
                chains[record['chain']].insert(len(chains[record['chain']])-1,line)

            if record['chain'] in [*chains.keys()]:
                chains[record['chain']].append(line)
                prev_res_no = record['res_no']
                prev_atom = record['atm_name']
            else:
                chains[record['chain']] = [line]
                prev_res_no= record['res_no']
                prev_atom = record['atm_name']



    return chains


def write_pdb(chainA, chainB, chain_break, outname):
    '''Save the CB coordinates for later processing
    '''

    with open(outname,'w') as file:
        #Write the first pdb file
        atomc=1
        aac=1
        chain = 'A'
        set_prev_res = True

        #Write chain A
        for line in chainA:
            record = parse_atm_record(line)

            if set_prev_res==True:
                prev_res_no=record['res_no']
                set_prev_res=False

            if record['res_no']>prev_res_no:
                aac+=1

            atom_blank = ' '*(5-len(str(atomc)))
            res_blank = ' '*(4-len(str(aac)))
            outline = line[:6]+atom_blank+str(atomc)+line[11:21]+chain+res_blank+str(aac)+line[26:]
            file.write(outline)
            atomc+=1

            prev_res_no = record['res_no']

        #Write chain B
        aac=1
        chain = 'B'
        set_prev_res = True
        for line in chainB:
            record = parse_atm_record(line)

            if set_prev_res==True:
                prev_res_no=record['res_no']
                set_prev_res=False

            if record['res_no']>prev_res_no:
                aac+=1

            atom_blank = ' '*(5-len(str(atomc)))
            res_blank = ' '*(4-len(str(aac)))
            outline = line[:6]+atom_blank+str(atomc)+line[11:21]+chain+res_blank+str(aac)+line[26:]
            file.write(outline)
            atomc+=1
            prev_res_no = record['res_no']


##################MAIN#######################

#Parse args
args = parser.parse_args()
#Data
pdbdir = args.pdbdir[0]
meta = pd.read_csv(args.meta[0])
outdir = args.outdir[0]

for i in range(len(meta)):
    print(i+1)
    row = meta.loc[i]
    pdbid = row.id1.split('_')[0]
    chains = read_all_chains(pdbdir+pdbid+'.pdb')
    chA, chB = row.id1.split('_')[1], row.id2.split('_')[1]
    #Write new chains
    write_pdb(chains[chA], chains[chB], row.l1, outdir+row.id1+'-'+row.id2+'.pdb')
