#!/bin/bash -x

### get script directory
BIN=$(dirname "${BASH_SOURCE[0]}")
pushd $BIN > /dev/null
BASEDIR=`pwd`
popd > /dev/null

### get input
PROT1=$1
PROT2=$2
CHAIN1=$3
CHAIN2=$4

### get prot IDs and define input format
SUFFIX=`echo $1 | awk -F . '{print $NF}'`
ID1=`basename ${PROT1%.$SUFFIX}`
ID2=`basename ${PROT2%.$SUFFIX}`

### singularity image
SINGULARITY=${BASEDIR}/AF_data/AF_environment.sif

### starting from structures
if [ $SUFFIX == "pdb" ] || [ $SUFFIX == "cif" ]; then
    ### extract proper chain (bakbone atoms only), rename it as A and renumber it from 1 not considering gaps
    singularity exec $SINGULARITY \
        python3 $BASEDIR/src/pdb_extract.py $PROT1 $CHAIN1 ./${ID1}_$CHAIN1.pdb
    singularity exec $SINGULARITY \
        python3 $BASEDIR/src/pdb_extract.py $PROT2 $CHAIN2 ./${ID2}_$CHAIN2.pdb

    ### rename chain2 as B
    sed 's/ A / B /g' ${ID2}_$CHAIN2.pdb > ${ID2}_${CHAIN2}_chainrename.pdb
    mv ${ID2}_${CHAIN2}_chainrename.pdb ${ID2}_$CHAIN2.pdb

    ### concatenate the two chains in a single structure for comparison with final model
    cat ${ID1}_$CHAIN1.pdb | grep -v 'END' > ${ID1}_${CHAIN1}-${ID2}_${CHAIN2}.pdb
    cat ${ID2}_$CHAIN2.pdb >> ${ID1}_${CHAIN1}-${ID2}_${CHAIN2}.pdb
    MERGEDSTR=${ID1}_${CHAIN1}-${ID2}_${CHAIN2}.pdb

    ### get single and concatenated fastas
    echo '>'${ID1}_${CHAIN1} > ${ID1}_$CHAIN1.fasta
    singularity exec $SINGULARITY \
        python3 $BASEDIR/src/sequence_extract.py ${ID1}_$CHAIN1.pdb A >> ${ID1}_$CHAIN1.fasta

    echo '>'${ID2}_${CHAIN2} > ${ID2}_$CHAIN2.fasta
    singularity exec $SINGULARITY \
	python3 $BASEDIR/src/sequence_extract.py ${ID2}_$CHAIN2.pdb B >> ${ID2}_$CHAIN2.fasta

    FASTA1=${ID1}_$CHAIN1.fasta
    FASTA2=${ID2}_$CHAIN2.fasta
    SEQ1=`tail -1 ${FASTA1}`
    SEQ2=`tail -1 ${FASTA2}`
    echo '>'${ID1}_${CHAIN1}-${ID2}_${CHAIN2} > ${ID1}_${CHAIN1}-${ID2}_${CHAIN2}.fasta
    echo ${SEQ1}${SEQ2} >> ${ID1}_${CHAIN1}-${ID2}_${CHAIN2}.fasta
    MERGED_FASTA=${ID1}_${CHAIN1}-${ID2}_${CHAIN2}.fasta
    MERGED_ID=${ID1}_${CHAIN1}-${ID2}_${CHAIN2}

### starting from sequences
elif [ $SUFFIX == "fasta" ]; then
    ### get merged sequences
    FASTA1=$PROT1
    FASTA2=$PROT2
    SEQ1=`grep -v '>' $PROT1 | tr -d '\n'`
    SEQ2=`grep -v '>' $PROT2 | tr -d '\n'`
    echo '>'${ID1}-${ID2} > ${ID1}-${ID2}.fasta
    echo ${SEQ1}${SEQ2} >> ${ID1}-${ID2}.fasta
    MERGED_FASTA=${ID1}-${ID2}.fasta
    MERGED_ID=${ID1}-${ID2}

### starting from msa
elif [ $SUFFIX == "a3m" ]; then
    head -2 $PROT1 > ${ID1}.fasta
    head -2 $PROT2 > ${ID2}.fasta
    FASTA1=${ID1}.fasta
    FASTA2=${ID2}.fasta
    SEQ1=`tail -1 $FASTA1`
    SEQ2=`tail -1 $FASTA2`
    echo '>'${ID1}-${ID2} > ${ID1}-${ID2}.fasta
    echo ${SEQ1}${SEQ2} >> ${ID1}-${ID2}.fasta
    MERGED_FASTA=${ID1}-${ID2}.fasta
    MERGED_ID=${ID1}-${ID2}
fi

### run alignment if not getting it as input
if [ $SUFFIX == "pdb" ] || [ $SUFFIX == "cif" ] || [ $SUFFIX == "fasta" ]; then
    singularity exec $SINGULARITY \
        hhblits -i $FASTA1 -d $BASEDIR/uniclust30_2018_08/uniclust30_2018_08 -E 0.001 -all -oa3m ${FASTA1%.fasta}.a3m
    A3M1=${FASTA1%.fasta}.a3m

    singularity exec $SINGULARITY \
        hhblits -i $FASTA2 -d $BASEDIR/uniclust30_2018_08/uniclust30_2018_08 -E 0.001 -all -oa3m ${FASTA2%.fasta}.a3m
    A3M2=${FASTA2%.fasta}.a3m

elif [ $SUFFIX == "a3m" ]; then
    A3M1=$PROT1
    A3M2=$PROT2
fi

### Unalign MSAs
singularity exec --nv --bind $BASEDIR:$BASEDIR $SINGULARITY \
    python3 $BASEDIR/src/unalign_MSA.py \
        --msa $A3M1 \
        --fasta unaligned_msa1.fasta

singularity exec --nv --bind $BASEDIR:$BASEDIR $SINGULARITY \
    python3 $BASEDIR/src/unalign_MSA.py \
        --msa $A3M2 \
        --fasta unaligned_msa2.fasta

### cluster unaligned fasta sequences at 62% sequence identity
$BASEDIR/src/cdhit/cd-hit -i unaligned_msa1.fasta -o clusters0.62_1 -c 0.62 -G 0 -n 3 -aS 0.9
$BASEDIR/src/cdhit/cd-hit -i unaligned_msa2.fasta -o clusters0.62_2 -c 0.62 -G 0 -n 3 -aS 0.9

### count clusters to obtain estimates of Neff scores
NEFF1=`grep '>Cluster' clusters0.62_1.clstr | wc -l`
NEFF2=`grep '>Cluster' clusters0.62_2.clstr | wc -l`

MGF=0.9
FUSEDMSA=${MERGED_ID}_fused.a3m
PAIREDMSA=${MERGED_ID}_paired.a3m
if [ -f $A3M1 ] && [ -f $A3M2  ]
then
    if [ ! -f $PAIREDMSA ]
    then
	python3 ${BASEDIR}/src/oxmatch.py --a3m1 $A3M1 --a3m2 $A3M2 --max_gap_fraction $MGF --outname $PAIREDMSA
    fi
    if [ ! -f $FUSEDMSA ]
    then
	python3 ${BASEDIR}/src/fuse_msas.py --a3m1 $A3M1 --a3m2 $A3M2 --max_gap_fraction $MGF --outname $FUSEDMSA
    fi
fi

### AlphaFold runtime parameters
L1=`${BASEDIR}/src/seqlen.sh $FASTA1` # Length of chain1 (where to introduce chain break)
L2=`${BASEDIR}/src/seqlen.sh $FASTA2` # Length of chain2
MSAS=${PAIREDMSA},${FUSEDMSA} # Comma separated list of msa paths
AFHOME=${BASEDIR}/src/alphafold/ # Path of alphafold directory in FoldDock
PARAM=${BASEDIR}/AF_data/  # Path to AF2 params
OUTFOLDER=./ # Path where AF2 generates its output folder structure
PRESET="full_dbs" # Configuration1 - no ensembling (full_dbs) and (reduced_dbs) or 8 model ensemblings (casp14).
MODEL_NAME="model_1" # Configuration2
MAX_RECYCLES=10 # recycles number (default=3)

# HEADER="if_plddt_av,if_plddt_std,ch1_plddt_av,ch1_plddt_std,ch2_plddt_av,ch2_plddt_std,"
# HEADER=$HEADER"plddt_av,plddt_std,num_atoms_in_interface,num_res_in_interface,Neff1,Neff2"
# if [ $SUFFIX == "pdb" ] || [ $SUFFIX == "cif" ]; then HEADER=$HEADER",DockQ"; fi
# echo $HEADER > ${MERGED_ID}_unranked.csv

for n in {1..5}; do
    cp $MERGED_FASTA ${MERGED_ID}_${n}.fasta

    ### Run Alphafold2 to fold provided chains
    singularity exec --nv --bind $BASEDIR:$BASEDIR $SINGULARITY \
        python3 $AFHOME/run_alphafold.py \
            --fasta_paths=${MERGED_ID}_${n}.fasta \
	    --output_dir=$OUTFOLDER \
	    --model_names=$MODEL_NAME \
            --max_recycles=$MAX_RECYCLES \
	    --data_dir=$PARAM \
	    --preset=$PRESET \
            --fold_only \
            --msas=$MSAS \
            --chain_break_list=$L1 \

    rm ${MERGED_ID}_${n}.fasta

  #   ### Run analysis on interface pLDDT and contacts
  #   INTSTATS=`singularity exec --nv --bind $BASEDIR:$BASEDIR $SINGULARITY \
  #                 python3 $BASEDIR/src/fetch_plDDT.py \
  #                     --model ${MERGED_ID}_$n/unrelaxed_model_1.pdb \
  #                     --metric ${MERGED_ID}_$n/result_model_1.pkl \
  #                     --l1 $L1 \
	#               --l2 $L2`
  #
  #   ### split modelled chains in 2, renumber and rename them A and B
  #   python3 $BASEDIR/src/split_chains.py --structure ${MERGED_ID}_$n/unrelaxed_model_1.pdb --outname ${MERGED_ID}_unranked_$n
  #   INTSTATS=${INTSTATS},${NEFF1},${NEFF2}
  #
  #   ### If starting from structures, compare with (merged) structures provided as input
  #   if [ $SUFFIX == "pdb" ] || [ $SUFFIX == "cif" ]; then
  #       DOCKQ=`singularity exec --nv --bind $BASEDIR:$BASEDIR $SINGULARITY \
	# 	   $BASEDIR/src/DockQ/DockQ.py \
	# 	       -short \
	# 	       ${MERGED_ID}_unranked_${n}.pdb \
	# 	       $MERGEDSTR | \
	# 	       grep "DockQ" | awk '{print $2}'`
  #
	# INTSTATS=${INTSTATS},$DOCKQ
  #   fi
  #
  #   echo $INTSTATS >> ${MERGED_ID}_unranked.csv
done

# ### rerank obtained structures
# if [ $SUFFIX == "pdb" ] || [ $SUFFIX == "cif" ]; then
#     python3 $BASEDIR/src/rerank.py --comm ${MERGED_ID}_unranked --csv ${MERGED_ID}_unranked.csv --col 9 --out ${MERGED_ID}_ranked
# else
#     python3 $BASEDIR/src/rerank.py --comm ${MERGED_ID}_unranked --csv ${MERGED_ID}_unranked.csv --col 9 --out ${MERGED_ID}_ranked
# fi
