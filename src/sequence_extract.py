import os
import sys
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M'}

def get_sequence(struc, chain, model=0):
    seq = ''
    for res in struc[model][chain].get_list():
        if not res.has_id("CA"): continue
        if res.get_resname() not in three2one:
            continue
        seq += three2one[res.get_resname()]

    return seq


if __name__ == '__main__':

    if not os.path.exists(sys.argv[1]):
        raise IOError ('Input structure {} not found!'.format(sys.argv[1]))
    if not sys.argv[1].endswith('.pdb') and not sys.argv[1].endswith('.cif'):
        raise IOError ('Structure must have a .pdb or .cif format.')

    pdbp = PDBParser(QUIET=True)
    cifp = MMCIFParser(QUIET=True)
    if sys.argv[1].endswith('.pdb'): struc = pdbp.get_structure('', sys.argv[1])
    else: struc = cifp.get_structure('', sys.argv[1])

    print (get_sequence(struc, sys.argv[2], model=0))
