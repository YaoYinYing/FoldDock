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
parser.add_argument('--atom_list', nargs=1, type= str, default=sys.stdin, help = 'List of atoms to parse.')
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

def read_all_chains_coords(pdbname, atom_list):
    '''Get all atom coordinates for all chains
    '''
    three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
    'SEC':'U', 'PYL':'O', 'GLX':'X', 'UNK': 'X'}
    with open(pdbname) as pdbfile:
        pdb_chains = {} #Coordinates
        sequence = {}
        prev_res_no=''
        prev_atm = ''
        for line in pdbfile:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            if record['atm_name']==prev_atm and record['res_no']==prev_res_no:
                continue
            if  record['atm_name'] not in atom_list:
                continue
            if record['chain'] in [*pdb_chains.keys()]:
                pdb_chains[record['chain']].append(line)
                prev_res_no= record['res_no']
                prev_atm = record['atm_name']
            else:
                pdb_chains[record['chain']] = [line]
                prev_res_no= record['res_no']
                prev_atm = record['atm_name']

            #Sequence
            if record['atm_name']=='CA':
                if record['chain'] in [*sequence.keys()]:
                        sequence[record['chain']]+=three_to_one[record['res_name']]
                else:
                    sequence[record['chain']] = three_to_one[record['res_name']]


    return pdb_chains, sequence

def write_fasta(id,sequence,outname):
    '''Write fasta files
    '''
    with open(outname, 'w') as file:
        file.write('>'+id+'\n')
        file.write(sequence)

def write_pdb(ch1, ch2, atom_list, outname):
    '''Save the extracted coordinates
    '''

    with open(outname,'w') as file:
        #Write the first pdb file
        atomc=1
        aac=1

        #Write chain 1
        chain = 'A'
        for line in ch1:
            record = parse_atm_record(line)
            atom_blank = ' '*(5-len(str(atomc)))
            res_blank = ' '*(4-len(str(aac)))
            outline = line[:6]+atom_blank+str(atomc)+line[11:21]+chain+res_blank+str(aac)+line[26:]
            file.write(outline)
            #Update index
            if atomc%len(atom_list)==0:
                aac+=1
            atomc+=1

        #Write chain 2
        chain = 'B'
        aac = 1
        for line in ch2:
            record = parse_atm_record(line)
            atom_blank = ' '*(5-len(str(atomc)))
            res_blank = ' '*(4-len(str(aac)))
            outline = line[:6]+atom_blank+str(atomc)+line[11:21]+chain+res_blank+str(aac)+line[26:]
            file.write(outline)
            #Update index
            if atomc%len(atom_list)==0:
                aac+=1
            atomc+=1




################MAIN###############
#Parse args
args = parser.parse_args()
#Data
pdbdir = args.pdbdir[0]
meta = pd.read_csv(args.meta[0])
atom_list = args.atom_list[0].split(',')
outdir = args.outdir[0]

#Extract all chains and write new pdbs, fasta and cat fasta
seqs1 = []
seqs2 = []
l1 = []
l2 = []
for i in range(len(meta)):
    row = meta.loc[i]
    pdbname = pdbdir+row.PDB+'.pdb'
    chains, sequences = read_all_chains_coords(pdbname, atom_list)
    try:
        ch1, seq1 = chains[row['Chain 1']], sequences[row['Chain 1']]
        ch2, seq2 = chains[row['Chain 2']], sequences[row['Chain 2']]
    except:
        pdb.set_trace()
    #Save
    seqs1.append(seq1)
    seqs2.append(seq2)
    l1.append(len(seq1))
    l2.append(len(seq2))

    #Write fasta
    write_fasta(row.PDB+'_'+row['Chain 1'],seq1,outdir+'fasta/'+row.PDB+'_'+row['Chain 1']+'.fasta')
    write_fasta(row.PDB+'_'+row['Chain 2'],seq2,outdir+'fasta/'+row.PDB+'_'+row['Chain 2']+'.fasta')
    write_fasta(row.PDB+'_'+row['Chain 1']+'-'+row.PDB+'_'+row['Chain 2'],seq1+seq2,outdir+'fasta/'+row.PDB+'_'+row['Chain 1']+'-'+row.PDB+'_'+row['Chain 2']+'.fasta')
    #Rewrite the files
    write_pdb(ch1, ch2, atom_list, pdbdir+row.PDB+'_'+row['Chain 1']+'-'+row.PDB+'_'+row['Chain 2']+'.pdb')


#Update df
meta['Seq1']=seqs1
meta['Seq2']=seqs2
meta['l1']=l1
meta['l2']=l2
meta.to_csv(args.meta[0])
