#Preprocess the new set
#Extract PDB
PDBDIR=./PDB/
META=./newset.csv
ATOM_LIST='N,CA,C'
OUTDIR=./
python3 ./extract_pdb.py --pdbdir $PDBDIR --meta $META --atom_list $ATOM_LIST --outdir $OUTDIR
