DSSPDIR=/home/pbryant/data/FoldDock/dssp/
CBDIR=./PDB/CB_coordinates/
META=./final_ids_lens.csv
OUTDIR=./
python3 ./parse_dssp.py --dsspdir $DSSPDIR --cbdir $CBDIR --meta $META --outdir $OUTDIR
