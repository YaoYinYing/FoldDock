import sys
import string
import argparse

parser = argparse.ArgumentParser(description = '''Split output structures resulting from AF2''')
parser.add_argument('--structure', type= str, default=sys.stdin, help = 'List of paths to msas to pair in a3m format.')
parser.add_argument('--outname', type= str, default='./', help = 'Path to file to write to.')
args = parser.parse_args()

sep = 0
chain = 'A'
prenum = 0
chains_ids = list(string.ascii_uppercase[1:])
with open(args.outname+'.pdb', 'w') as out:
    for line in open(str(args.structure)):
        if not line.startswith('ATOM'): continue
        resnum = int(line[22:26].strip())
        if resnum > prenum+1: 
            sep = resnum-1
            chain = chains_ids.pop(0)
        
        newnum = resnum - sep
        out.write(line[:21]+chain+str(newnum).rjust(4)+line[26:])
        prenum = resnum

