import argparse
import sys
import os
import pandas as pd
import pdb

parser = argparse.ArgumentParser(description = '''Concatenate the sequences of two fasta files into 1.''')

parser.add_argument('--fastadir', nargs=1, type= str, default=sys.stdin, help = 'Path to data.')
parser.add_argument('--meta', nargs=1, type= str, default=sys.stdin, help = 'Path to csv with information on what sequences to merge.')
parser.add_argument('--mode', nargs=1, type= str, default=sys.stdin, help = 'Mode: bench4 or marks.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output directory. Include /in end')


def read_fasta(infile):
    '''Read fasta'''
    sequence=''
    with open(infile, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('>'): #OX=OrganismIdentifier
                continue
            else:
                sequence += line


    return sequence

def write_fasta(meta, mode, outdir):
    '''Write the merged fasta sequences
    '''
    for i in range(len(meta)):
        row = meta.loc[i]

        if mode=='bench4':
            with open(outdir+row.PDB+'_u1-'+row.PDB+'_u2.fasta','w') as file:
                file.write('>'+row.PDB+'_u1-'+row.PDB+'_u2\n')
                file.write(row.concatenated_seq)
        if mode=='marks':
            with open(outdir+row.id1+'-'+row.id2+'.fasta','w') as file:
                file.write('>'+row.id1+'-'+row.id2+'\n')
                file.write(row.concatenated_seq)
        if mode=='CASP14':
            id1, id2 = row.PDB.split('-')
            with open(outdir+id1+'-'+id2+'.fasta','w') as file:
                file.write('>'+id1+'-'+id2+'\n')
                file.write(row.concatenated_seq)
        if mode=='marks_neg':
            id1, id2 = row['Uniprot ID 1'], row['Uniprot ID 2']
            with open(outdir+id1+'-'+id2+'.fasta','w') as file:
                file.write('>'+id1+'-'+id2+'\n')
                file.write(row.concatenated_seq)
        if mode=='negatome':
            with open(outdir+row.id1+'-'+row.id2+'.fasta','w') as file:
                file.write('>'+row.id1+'-'+row.id2+'\n')
                file.write(row.concatenated_seq)
#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
fastadir = args.fastadir[0]
meta = pd.read_csv(args.meta[0])
mode = args.mode[0]
outdir = args.outdir[0]

#Results
sequences_1 = []
sequences_2 = []
merged_sequences = []
l1 = []
l2 = []
for i in range(len(meta)):
    row = meta.loc[i]
    if mode=='bench4':
        seq1 = read_fasta(fastadir+row.PDB+'_u1.fasta')
        seq2 = read_fasta(fastadir+row.PDB+'_u2.fasta')
    if mode=='marks':
        seq1 = read_fasta(fastadir+row.id1+'.fasta')
        seq2 = read_fasta(fastadir+row.id2+'.fasta')
    if mode=='CASP14':
        id1, id2 = row.PDB.split('-')
        seq1 = read_fasta(fastadir+id1+'.fasta')
        seq2 = read_fasta(fastadir+id2+'.fasta')
    if mode=='marks_neg':
        id1, id2 = row['Uniprot ID 1'], row['Uniprot ID 2']
        seq1 = read_fasta(fastadir+id1+'.fasta')
        seq2 = read_fasta(fastadir+id2+'.fasta')
    if mode=='negatome':
        seq1 = read_fasta(fastadir+row.id1+'.fasta')
        seq2 = read_fasta(fastadir+row.id2+'.fasta')

    #Save
    sequences_1.append(seq1)
    sequences_2.append(seq2)
    merged_sequences.append(seq1+seq2)
    l1.append(len(seq1))
    l2.append(len(seq2))

#Add to meta
meta['seq1']=sequences_1
meta['seq2']=sequences_2
meta['concatenated_seq']=merged_sequences
meta['l1']=l1
meta['l2']=l2
meta.to_csv(args.meta[0])

#Write merged
write_fasta(meta, mode, outdir)
