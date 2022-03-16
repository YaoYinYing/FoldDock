import sys
import numpy as np
import subprocess as sp
import argparse

parser = argparse.ArgumentParser(description = '''Split output structures resulting from AF2''')
parser.add_argument('--comm', type= str, default=sys.stdin, help = 'Common path for 5 AF2 models with suffix _N.pdb (with N being an int from 1 to 5).')
parser.add_argument('--csv', type= str, default=sys.stdin, help = 'Path to csv file with DockQ and pLDDT scores')
parser.add_argument('--col', type= int, default=sys.stdin, help = 'csv column to use for model ranking')
parser.add_argument('--out', type= str, default=sys.stdin, help = 'Common path + prefix for output ranked models')
args = parser.parse_args()

with open(args.csv) as csv: header = csv.readline()

scoremat = [[float(val) for val in line.split(',')] for line in open(args.csv) \
            if not line.startswith('if_plddt_av')]
scoremat = np.array(scoremat)
order = np.argsort(scoremat[:,args.col-1])[::-1]

with open(args.out+'.csv','w') as out:
    out.write('rank,'+header)
    for idx, n in enumerate(order):
        sp.run('cp {} {}'.format(args.comm+'_'+str(n+1)+'.pdb', args.out+'_'+str(idx+1)+'.pdb') ,shell=True)
        line = ','.join([str(val) for val in scoremat[n]])
        out.write(str(idx+1)+','+line+'\n')
