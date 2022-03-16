BASE=/proj/berzelius-2021-29/users/x_patbr #Path to home
REWRITE=$BASE/FoldDock/src/analysis/dockq/rewrite_af_pdb.py
##Rewrite the AF2 output to contain 2 chains instead of 1 for the DockQ evaluation

################################
###########new_dimers###########
################################

##########AF std and HHblits n2 msas############
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/new_dimers/af_std_and_hhblits_msas/model_1/recycle_10
META=$BASE/FoldDock/data/new_dimers/newset.csv
MODE='new_dimers'

for i in {1..5}
do
  OUTDIR=$PDBDIR/run_$i/rewritten_pdbs_all_atoms/
  mkdir $OUTDIR
  python3 $REWRITE --pdbdir $PDBDIR/run_$i/ --meta $META --mode $MODE  --outdir $OUTDIR
done


################################
############Marks###############
################################
META=$BASE/FoldDock/data/marks/final_ids_lens.csv
MODE='marks'
##########HHblits n2##########
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/marks/hhblits_msa/model_1/recycle_10/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR


##########AF std msa############
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/marks/af2_std_msa/model_1/recycle_10/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR

##########HHblits n2: paired+fused##########
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/marks/fused_and_hhbllits_msas/model_1/recycle_10
for i in {1..5}
do
  OUTDIR=$PDBDIR'/run_'$i'/rewritten_pdbs_all_atoms/'
  mkdir $OUTDIR
  python3 $REWRITE --pdbdir $PDBDIR/run_$i/ --meta $META --mode $MODE  --outdir $OUTDIR
done

##########AF std and HHblits n2 msas############
#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/marks/af_std_and_hhblits_msas/model_1/recycle_10

for i in {1..5}
  do
  OUTDIR=$PDBDIR'/run_'$i'/rewritten_pdbs_all_atoms/'
  mkdir $OUTDIR
  python3 $REWRITE --pdbdir $PDBDIR/run_$i/ --meta $META --mode $MODE  --outdir $OUTDIR
done


#############################
########Benchmark 4##########
#############################
META=$BASE/FoldDock/data/dockground/dockground.csv
MODE='bench4'
##########HHblits n2##########
#model_1_ptm
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/hhblits_msa/model_1_ptm/recycle_10/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR


#model_1_ptm
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/hhblits_msa/model_1_ptm/ensemble_8/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR

#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/hhblits_msa/model_1/recycle_10/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR


#model_1
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/hhblits_msa/model_1/ensemble_8/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR


##########AF std MSA##########
#model_1_ptm
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af2_std_msa/model_1_ptm/recycle_10/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR

#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af2_std_msa/model_1_ptm/ensemble_8/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR


#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af2_std_msa/model_1/recycle_10/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR

#model_1
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af2_std_msa/model_1/ensemble_8/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR


##########AF std and HHblits n2 msas############
#model_1_ptm
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af_std_and_hhblits_msas/model_1_ptm/recycle_10/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR

#model_1_ptm
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af_std_and_hhblits_msas/model_1_ptm/ensemble_8/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR

#model_1
#reycle 10
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af_std_and_hhblits_msas/model_1/recycle_10/run_1/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR


#model_1
#ensemble 8
PDBDIR=/proj/berzelius-2021-29/users/x_patbr/results/FoldDock/bench4/af_std_and_hhblits_msas/model_1/ensemble_8/
OUTDIR=$PDBDIR/rewritten_pdbs_all_atoms/
mkdir $OUTDIR
python3 $REWRITE --pdbdir $PDBDIR --meta $META --mode $MODE  --outdir $OUTDIR
