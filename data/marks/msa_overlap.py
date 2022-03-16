import argparse
import sys
import os
import numpy as np
from collections import Counter
import pandas as pd
import re
import pdb

parser = argparse.ArgumentParser(description = '''Match OX identifiers for interacting chains.''')

parser.add_argument('--filename', nargs=1, type= str, default=sys.stdin, help = 'Path to a3m file.')
parser.add_argument('--l1', nargs=1, type= int, default=sys.stdin, help = 'Length of chain 1.')
parser.add_argument('--l2', nargs=1, type= int, default=sys.stdin, help = 'Length of chain 2.')
parser.add_argument('--outname', nargs=1, type= str, default=sys.stdin, help = 'Outname')



def get_overlap(filename, l1, l2):
    '''Analyze the overlap between chains from a cat chain search
    '''
    overlap = []
    with open(filename,'r') as file:
        for line in file:
            if line[0]=='>':
                continue
            else:
                line = line.rstrip()
                #Read only upper case chars - lower case are inserts
                line = ''.join(re.findall('[A-Z,-]',line))
                if len(line)!=l1+l2:
                    continue
                try:
                    g1 = Counter(line[:l1])['-']/l1
                except:
                    g1=0
                try:
                    g2 = Counter(line[l1:])['-']/l2
                except:
                    g2=0

                #Save line
                if g1==g2:
                    continue
                if g1>g2: #More gaps in beginning --> bottom
                    #See how many % are non gaps in g1
                    overlap.append(1-g1)
                else:
                    overlap.append(1-g2)

    return np.average(overlap)

#################MAIN####################
#Parse args
args = parser.parse_args()
#Data
filename = args.filename[0]
l1 = args.l1[0]
l2 = args.l2[0]
outname = args.outname[0]
overlap = get_overlap(filename, l1, l2)
df = pd.DataFrame()
df['file']=[filename]
df['overlap']=[overlap]
df.to_csv(outname, header=None)
