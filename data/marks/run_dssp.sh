#Run DSSP to determine the 2ndary structure of protein complexes

DSSP=/home/pbryant/dssp
for file in /home/pbryant/data/MSAmatch/pdb/*.pdb
do
	echo $file
	$DSSP -i $file -o $file.dssp
done
