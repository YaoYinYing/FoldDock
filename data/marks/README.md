## Data
A set of 1,661 protein complexes with known interfaces from https://www.nature.com/articles/s41467-021-21636-z has been used here. These complexes share less than 30 % sequence identity, have a resolution between 1-5 Å and constitute unique pairs of PFAM domains (no single protein pair have PFAM domains matching that of any other pair).

## PDB files
To download these complexes do: cd ./PDB
bash ./batch_download.sh -f full_set_unique_ids.txt -p

Then extract all fasta sequences and CB coordinates found in these PDB files accordingly:

cd PDB
PDBDIR=#path to pdb directory with downloaded files
META=../full_set.csv
OUTDIR=../fasta/
python3 ./pdb_parser.py --pdbdir $PDBDIR --pdbmeta $META --outdir $OUTDIR

Note that only 1554/1661 out of these complexes contain CBs.

## HHblits
The fasta files are then used to run hhblits (version 3.1.0) against uniclust30_2018_08 using the criteria "-E 0.001 -all -oa3m"

From the resulting .a3m files, different matching strategies are compared by comparing selections of orthologs based on the "OX" (organism) identifiers in the a3m files. The strategies compared are 1. taking the "top" hits, 2. taking the "bottom" hits and 3. simply concatenating the .a3m files (no matching).

All of these steps are outlined in ./hhblits/preprocess.sh (and match_ox.py). The number of possible matches, number of mathing OX and more were analyzed in ./eda (./eda/eda.sh).

## Protein features
Several protein features have been extracted from the complexes:

- Alignment scores between the chains in each complex

- Secondary structure annotations from DSSP
