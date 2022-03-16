#Fetch the atos N, CA and C to do the dockq comparisons
INDIR=./PDB
OUTDIR=./PDB/rewritten_complex/
for i in {1..219}
do
  ID=$(sed -n $i'p' pdb_ids.txt)
  FILE=$INDIR/$ID'_bc.pdb'
  echo $i
  grep 'N   \|CA \|C  ' $FILE > $OUTDIR/$ID'_bc.pdb'

done
