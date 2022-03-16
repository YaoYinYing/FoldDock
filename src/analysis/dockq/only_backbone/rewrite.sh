

##Rewrite the AF2 output to contain 2 chains instead of 1 for the DockQ evaluation

################################
############new_dimers###############
################################

##########AF std and HHblits n2 msas############
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/new_dimers/af_std_and_hhblits_msas/model_1/recycle_10
META=../data/new_dimers/newset.csv
MODE='new_dimers'
ATOMS='N,CA,C'
#for i in {1..5}
#do
  #OUTDIR=$PDBDIR/run_$i/rewritten_pdbs/
  #mkdir $OUTDIR
  #python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR/run_$i/ --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR
#done


################################
############Marks###############
################################

##########HHblits n2##########
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/hhblits_msa/model_1/recycle_10/
META=../data/marks/final_ids_lens.csv
MODE='marks'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


##########AF std msa############
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/af2_std_msa/model_1/recycle_10/
META=../data/marks/final_ids_lens.csv
MODE='marks'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR

##########HHblits n2: paired+fused##########
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/fused_and_hhbllits_msas/model_1/recycle_10
META=../data/marks/final_ids_lens.csv
MODE='marks'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR
#for i in {2..5}
#do
#mkdir $OUTDIR/run_$i/rewritten_pdbs/
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR/run_$i/ --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR/run_$i/rewritten_pdbs/
#done

##########AF std and HHblits n2 msas############
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/af_std_and_hhblits_msas/model_1/recycle_10
META=../data/marks/final_ids_lens.csv
MODE='marks'
ATOMS='N,CA,C'
# for i in {2..5}
# do
# OUTDIR=$PDBDIR/run_$i/rewritten_pdbs/
# mkdir $OUTDIR
# python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR/run_$i/ --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR
# done


#############################
########Benchmark 4##########
#############################

##########HHblits n2##########
#model_1_ptm
#reycle 10
PDBDIR=/proj/snic2019-35-62/users/x_patbr/results/af2/hhblits_n2_merged/model_1_ptm/recycle_10/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


#model_1_ptm
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/hhblits_msa/model_1_ptm/ensemble_8/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
# OUTDIR=$PDBDIR/rewritten_pdbs/
# mkdir $OUTDIR
# python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


#model_1
#reycle 10
#PDBDIR=/proj/snic2019-35-62/users/x_patbr/results/af2/hhblits_n2_merged/model_1/recycle_10/
#META=../data/dockground/dockground.csv
MODE='bench4'
#ATOMS='N,CA,C'
#OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


#model_1
#ensemble 8
PDBDIR=/proj/snic2019-35-62/users/x_patbr/results/af2/hhblits_n2_merged/model_1/ensemble_8/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


##########AF std MSA##########
#model_1_ptm
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af2_std_msa/model_1_ptm/recycle_10/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
# mkdir $OUTDIR
# python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR

#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af2_std_msa/model_1_ptm/ensemble_8/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af2_std_msa/model_1/recycle_10/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR

#model_1
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af2_std_msa/model_1/ensemble_8/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


#model_2
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af2_std_msa/model_2/ensemble_8/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
mkdir $OUTDIR
python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR


##########AF std and HHblits n2 msas############
#model_1_ptm
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af_std_and_hhblits_msas/model_1_ptm/recycle_10/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR

#model_1_ptm
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af_std_and_hhblits_msas/model_1_ptm/ensemble_8/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR

#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af_std_and_hhblits_msas/model_1/recycle_10/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
#OUTDIR=$PDBDIR/rewritten_pdbs/
#mkdir $OUTDIR
#python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR

#Run 1-5
# for i in {1..5}
# do
# 	OUTDIR=$PDBDIR'run_'$i'/rewritten_pdbs/'
# 	mkdir $OUTDIR
# 	python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR'run_'$i'/' --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR
# done


#model_1
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af_std_and_hhblits_msas/model_1/ensemble_8/
META=../data/dockground/dockground.csv
MODE='bench4'
ATOMS='N,CA,C'
OUTDIR=$PDBDIR/rewritten_pdbs/
# mkdir $OUTDIR
# python3 ./rewrite_af_pdb.py --pdbdir $PDBDIR --meta $META --mode $MODE --atom_list $ATOMS --outdir $OUTDIR
