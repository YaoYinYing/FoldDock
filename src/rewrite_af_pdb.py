import argparse
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Rewrite the AF2 output to contain two different chains in the PDB files.''')

parser.add_argument('--pdbfile', nargs=1, type= str, default=sys.stdin, help = 'Path to predicted pdbfile.')
parser.add_argument('--l1', nargs=1, type= int, default=sys.stdin, help = 'Length of first chain in pdbfile.')
parser.add_argument('--outname', nargs=1, type= str, default=sys.stdin, help = 'Path to rewritten pdbfile')

##########Functions##############
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

        for line in pdbfile:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            if record['chain'] in [*pdb_chains.keys()]:
                pdb_chains[record['chain']].append(line)
            else:
                pdb_chains[record['chain']] = [line]


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
pdbfile = args.pdbfile[0]
l1 = args.l1[0]
outname = args.outname[0]

#Read PDB
chains = read_all_chains_coords(pdbfile)
#Rewrite the files
write_pdb(chains['A'], l1, outname)
