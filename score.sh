#This is a simple scoring script to calculate pDockQ
#1 Rewrite the predicted pdbfile to contain two chains (A and B)
PDBDIR=./test/ #Path to predicted structure
ID=1ay7 #Name of predicted structure
L1=96 #Length of first chain of predicted structure
#Rewrite
python3 ./src/rewrite_af_pdb.py --pdbfile $PDBDIR/unrelaxed_model_1.pdb --l1 $L1 --outname $PDBDIR/$ID'_pred.pdb'

#2. Calculate the predicted DockQ score (pDockQ)
#pDockQ
python3 ./src/pdockq.py --pdbfile $PDBDIR/$ID'_pred.pdb' --pickle_file $PDBDIR/result_model_1.pkl
