import sys
import glob
import argparse
import numpy as np
from Bio import AlignIO

parser = argparse.ArgumentParser(description = '''Split output structures resulting from AF2''')
parser.add_argument('--msa', type= str, default=sys.stdin, help = 'Input MSA')
parser.add_argument('--fasta', type= str, default=sys.stdin, help = 'Output unaligned FASTA file')
args = parser.parse_args()

with open(args.fasta, 'w') as out:
    if args.msa[-3:] == 'sto':
        alignment = AlignIO.read(args.msa, 'stockholm')
        fasta_alignment = format(alignment, 'fasta')
        
    elif args.msa[-3:] == 'a3m':
        fasta_alignment = ''.join([line for line in open(args.msa)])
    
    else: raise TypeError('Only Stockholm and A3M formats are supported!')
    
    seq = ''
    for line in fasta_alignment.split('\n'):
        if line.startswith('>'): 
            if seq != '': out.write(key+'\n'+seq+'\n')
            key = line.strip()
            seq = ''

        else: seq += ''.join([char for char in line.strip()\
                             if char != '-' and char.isupper()])
