import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Read an a3m file and get all "OX" identifiers.''')

parser.add_argument('--a3mdir', nargs=1, type= str, default=sys.stdin, help = 'Path to data.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')


def read_a3m(infile,max_gap_fraction=0.9):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
             'G': 6,'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,'N': 12,
             'O': 21, 'P': 13,'Q': 14, 'R': 15, 'S': 16, 'T': 17,
             'V': 18, 'W': 19, 'Y': 20,'U': 21, 'Z': 21, 'X': 21, 'J': 21}

    parsed = []#Save extracted msa
    species = []
    seqlen = 0
    lc = 0
    with open(infile, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('>'): #OX=OrganismIdentifier
                if 'OX=' in line:
                    OX= line.split('OX=')[1]
                    if len(OX)>0:
                        species.append(int(OX.split(' ')[0]))
                    else:
                        species.append(0)
                else:
                    species.append(0)
                continue
            line = line.rstrip()
            gap_fraction = line.count('-') / float(len(line))
            if gap_fraction <= max_gap_fraction:#Only use the lines with less than 90 % gaps
                parsed.append([mapping.get(ch, 22) for ch in line if not ch.islower()])
            else:
                if len(species)>1:
                    species = species[:-1] #Remove the previously stored species
                    continue
            #Check that the lengths match
            if len(parsed[-1])!=seqlen and lc>=1:
                parsed = parsed[:-1]
                species = species[:-1]
                continue
            seqlen = len(parsed[-1])
            lc+=1


    return np.array(parsed, dtype=np.int8, order='F'), np.array(species)


#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
a3mdir = args.a3mdir[0]
outdir = args.outdir[0]

a3mfiles = glob.glob(a3mdir+'*.a3m')
#Results
ids = []
length = []
nseqs_total = []
nseqs_OX = []
nunique_OX = []
for name in a3mfiles:
    ids.append(name.split('/')[-1].split('.')[0])
    msa, species = read_a3m(name)
    length.append(msa.shape[1])
    nseqs_total.append(msa.shape[0])
    nseqs_OX.append(np.argwhere(species!=0).shape[0])
    nunique_OX.append(np.unique(species).shape[0]-1) #Subtract 1 (no species)
    #Save the species to compute stats
    np.save(outdir+ids[-1]+'_species.npy',species)
    #Save the a3m matrix
    np.save(outdir+ids[-1]+'_a3m_matrix.npy',msa)



#DF
results = pd.DataFrame()
results['ID'] = ids
results['Seqlen'] = length
results['nseqs_total'] = nseqs_total
results['nseqs_with_OX'] = nseqs_OX
results['nunique_OX'] = nunique_OX
results['OX_fraction_seqs'] = results.nseqs_with_OX/results.nseqs_total
results['OX_fraction_unique'] = results.nunique_OX/results.nseqs_with_OX
results.to_csv(outdir+'a3mstats.csv', index=None)
