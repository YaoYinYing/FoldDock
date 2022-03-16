import os
import sys
import glob
import string
from typing import Optional

from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import Select
from Bio.PDB.PDBIO import PDBIO

select_atoms=['N', 'CA', 'C']

class Select_Atoms(Select):
    def accept_atom(self, atom):
        return atom.get_id() in select_atoms

def extract(struc: Structure, 
            chains: list,
            outpath: Optional[str],
            select_atoms: Optional[list],
            renumber: Optional[bool],
            rename: Optional[bool],
            fill: Optional[bool]) -> Optional[Structure]:

    if outpath:
        if not os.path.exists(os.path.dirname(outpath)):
            raise IOError ('Outout destination {} not found!'.format(outpath))
        if not outpath.endswith('.pdb'):
            raise IOError ('Specify full path with .pdb format for output structure.')

    for chain in chains:
        if len(chain) != 2 and len(chain) != 4:
            raise IOError ('Each element in chains must be either a \
                            2 elements [model_id, chain_id] or 4 elements \
                            lists [model_id, chain_id, first_residue_id, \
                            last_residue_id]. Faulty element:{}'.format(chain))
        if (type(chain[0]) != int):
            raise IOError ('First element in each list in chains must be a \
                            integer indicating model_id. Faulty element:{}'.format(chain))
        if (type(chain[1]) != str):
            raise IOError ('Second element in each list in chains must be a \
                            string indicating chain_id. Faulty element:{}'.format(chain))
        if len(chain) == 4:
            if type(chain[2]) != int or type(chain[3]) != int:
                raise IOError ('Third and fourth elements in each list in chains \
                                must be integers indicating first and last residues \
                                to extract. Faulty element:{}'.format(chain))

    available_ids = list(string.ascii_uppercase)+list(string.ascii_lowercase)

    iopdb = PDBIO()
    custom_model = Model(0)
    custom_struc = Structure('xxxx')
    
    if chains == []:
        chains = [[model.get_id(), chain.get_id()] \
                  for model in struc for chain in model]
        
    print (chains)
    for target in chains:
        model = target[0]
        chain = target[1]
        
        target_chain = struc[model][chain]
        target_chain.detach_parent()

        if renumber:
            residues = unfold_entities(target_chain, 'R')
            target_chain = Chain(target_chain.get_id())
            delta = residues[0].get_id()[1]-1
            for residue in residues:
                if residue.get_id()[0] != ' ': continue
                residue.detach_parent()
                new = residue.get_id()[1] - delta
                residue.id = (' ', new, ' ')
                target_chain.add(residue)

        if fill:
            residues = unfold_entities(target_chain, 'R')
            for idx, residue in enumerate(residues):
                residue.id = (' ', idx+1, ' ')

        if rename:
            target_chain.id = available_ids.pop(0)
            custom_model.add(target_chain)

        if len(target) > 2:
            raise NotImplementedError

    custom_struc.add(custom_model)
    iopdb.set_structure(custom_struc)

    if outpath:
        if select_atoms: iopdb.save(outpath, Select_Atoms())
        else: iopdb.save(outpath)
    else: 
        if select_atoms: raise NotImplementedError
        else: return custom_struc


def main():
    chainlist = [[0, sys.argv[2].split('_')[0]]]

    if not os.path.exists(sys.argv[1]):
        raise IOError ('Input structure {} not found!'.format(sys.argv[1]))
    if not sys.argv[1].endswith('.pdb') and not sys.argv[1].endswith('.cif'):
        raise IOError ('Structure must have a .pdb or .cif format.')

    pdbp = PDBParser(QUIET=True)
    cifp = MMCIFParser(QUIET=True)
    if sys.argv[1].endswith('.pdb'): struc = pdbp.get_structure('', sys.argv[1])
    else: struc = cifp.get_structure('', sys.argv[1])

    extract(struc, 
            chainlist,
            sys.argv[3],
            select_atoms=['N', 'CA', 'C'],
            renumber=True,
            rename=True,
            fill=True)


if __name__ == '__main__':
    main()

