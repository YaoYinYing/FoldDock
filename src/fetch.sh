#!/bin/bash -x
#Fetch plDDT from predicted structures and analyze results

dir=/proj/berzelius-2021-29/users/x_arnel/Hu.Map-2.0/
#########New dimers###########
DIR=$dir/pdb/$1-$2/
MODEL=$DIR/unrelaxed_model_1.pdb
METRIC=$DIR/result_model_1.pkl
META=../data/new_dimers/newset.csv
IT=10 #Interface threshold Å
FETCH_ATOMS='CB'
CBR=200 #How many residues that were used as a chain break insertion

L1=`$dir/bin/seqlen.bash $dir/seq/$1.fasta`
L2=`$dir/bin/seqlen.bash $dir/seq/$2.fasta`
#for i in {1..5}
#do
python3 $dir/bin/fetch_plDDT.py --model $MODEL --metric $METRIC --it $IT --fetch_atoms $FETCH_ATOMS --cbr $CBR --l1 $L1 --l2 $L2
#done

