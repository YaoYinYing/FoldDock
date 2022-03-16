
#First run the a3m_parser for all your a3m files to obtain a list of
#organism identifiers and an array of the MSA
#Parse a3m
A3MDIR=#Path to a3m files
OUTDIR=./
#python3 ./a3m_parser.py --a3mdir $A3MDIR --outdir $OUTDIR


#Match OX (organism identifiers) and wrte statistics
OXDIR=../data/hhblits/
META=../data/full_set.csv
OUTDIR=../data/hhblits/
#python3 ./match_ox.py --oxdir $OXDIR --pdbmeta $META --outdir $OUTDIR

#Select sequences based on OX
OXDIR=#Path to directory with species information and a3m matrices in np format
META=../full_set_kingdom.csv
MODE=top #bottom or all = no selection
IDLIST=#Optionally procid a list of ids, otherwise all in meta will be used
OUTDIR=./
#python3 ./select_ox.py --oxdir $OXDIR --pdbmeta $META --mode $MODE --idlist $IDLIST --outdir $OUTDIR


#Calculate the Neff using an 80% similarity cutoff
A3MDIR=/home/patrick/MSAMatch/data/hhblits/testfiles/ #Path to a3m files
META=../full_set_kingdom.csv
MODE='top'
OUTDIR=./
python3 ./neff.py --a3mdir $A3MDIR --pdbmeta $META --mode $MODE --outdir $OUTDIR
