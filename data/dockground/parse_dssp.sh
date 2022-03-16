DSSPDIR=./PDB/DSSP/
PDBDIR=./PDB/
META=./dockground.csv
OUTDIR=./
python3 ./parse_dssp.py --dsspdir $DSSPDIR --pdbdir $PDBDIR --meta $META --outdir $OUTDIR
