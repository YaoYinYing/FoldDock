#Fetch plDDT from predicted structures and analyze results

#########New dimers###########
MODELDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/new_dimers/af_std_and_hhblits_msas/model_1/recycle_10/
METRICDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/new_dimers/af_std_and_hhblits_msas/model_1/recycle_10/
META=../data/new_dimers/newset.csv
IT=8 #Interface threshold Å
FETCH_ATOMS='CB'
CBR=200 #How many residues that were used as a chain break insertion
DFSUFFIX='newset_af_std_and_hhblits_msas' #Suffix for results df
MODE='newset'
OUTDIR=./
#for i in {1..5}
#do
#python3 ./fetch_plDDT.py --modeldir $MODELDIR'/run_'$i'/' --metricdir $METRICDIR'/run_'$i'/' --meta $META --it $IT --fetch_atoms $FETCH_ATOMS --cbr $CBR --df_suffix $DFSUFFIX$i --mode $MODE --outdir $OUTDIR
#done

##########Marks set##########
#AF2+Paired MSAs
MODELDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/af_std_and_hhblits_msas/model_1/recycle_10/
METRICDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/af_std_and_hhblits_msas/model_1/recycle_10/
META=../data/marks/final_ids_lens.csv
IT=8 #Interface threshold Å
FETCH_ATOMS='CB'
CBR=200 #How many residues that were used as a chain break insertion
DFSUFFIX='marks_af_std_and_hhblits_msas' #Suffix for results df
MODE='marks'
OUTDIR=./
# for i in {2..5}
# do
# python3 ./fetch_plDDT.py --modeldir $MODELDIR'/run_'$i'/' --metricdir $METRICDIR'/run_'$i'/' --meta $META --it $IT --fetch_atoms $FETCH_ATOMS --cbr $CBR --df_suffix $DFSUFFIX$i --mode $MODE --outdir $OUTDIR
# done


#Fused+Paired MSAs
MODELDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/fused_and_hhbllits_msas/model_1/recycle_10/
METRICDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks/fused_and_hhbllits_msas/model_1/recycle_10/
META=../data/marks/final_ids_lens.csv
IT=8 #Interface threshold Å
FETCH_ATOMS='CB'
CBR=200 #How many residues that were used as a chain break insertion
DFSUFFIX='marks_fused_and_hhblits_msas' #Suffix for results df
MODE='marks'
OUTDIR=./
#for i in {1..5}
#do
  #python3 ./fetch_plDDT.py --modeldir $MODELDIR'/run_'$i'/' --metricdir $METRICDIR'/run_'$i'/' --meta $META --it $IT --fetch_atoms $FETCH_ATOMS --cbr $CBR --df_suffix $DFSUFFIX$i --mode $MODE --outdir $OUTDIR
#done


##########Marks negative set##########
#Fused+Paired MSAs
MODELDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks_negative/fused_and_hhbllits_msas/model_1/recycle_10/
METRICDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/marks_negative/fused_and_hhbllits_msas/model_1/recycle_10/
META=../data/marks/marks_negative_set/neg_set_ids.csv
IT=8 #Interface threshold Å
FETCH_ATOMS='CB'
CBR=200 #How many residues that were used as a chain break insertion
DFSUFFIX='marks_negative_fused_and_hhblits_msas' #Suffix for results df
MODE='marks'
OUTDIR=./
#python3 ./fetch_plDDT.py --modeldir $MODELDIR --metricdir $METRICDIR --meta $META --it $IT --fetch_atoms $FETCH_ATOMS --cbr $CBR --df_suffix $DFSUFFIX$i --mode $MODE --outdir $OUTDIR

##########Negatome##########
#Fused+Paired MSAs
MODELDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/negatome/fused_and_hhbllits_msas/model_1/recycle_10/
METRICDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/negatome/fused_and_hhbllits_msas/model_1/recycle_10/
META=../data/negatome/negatome_ids.csv
IT=8 #Interface threshold Å
FETCH_ATOMS='CB'
CBR=200 #How many residues that were used as a chain break insertion
DFSUFFIX='negatome_fused_and_hhblits_msas' #Suffix for results df
MODE='marks'
OUTDIR=./
python3 ./fetch_plDDT.py --modeldir $MODELDIR --metricdir $METRICDIR --meta $META --it $IT --fetch_atoms $FETCH_ATOMS --cbr $CBR --df_suffix $DFSUFFIX$i --mode $MODE --outdir $OUTDIR



#########Benchmark 4###########
MODELDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af_std_and_hhblits_msas/model_1/recycle_10/
METRICDIR=/proj/berzelius-2021-29/users/x_patbr/results/af2/af_std_and_hhblits_msas/model_1/recycle_10/
META=../data/dockground/dockground.csv
IT=8 #Interface threshold Å
FETCH_ATOMS='CB'
CBR=200 #How many residues that were used as a chain break insertion
DFSUFFIX='bench4_af_std_and_hhblits_msas_model_1_recycle_10_run_' #Suffix for results df
MODE='bench4'
OUTDIR=./
#for i in {1..5}
#do
#python3 ./fetch_plDDT.py --modeldir $MODELDIR'/run_'$i'/' --metricdir $METRICDIR'/run_'$i'/' --meta $META --it $IT --fetch_atoms $FETCH_ATOMS --cbr $CBR --df_suffix $DFSUFFIX$i --mode $MODE --outdir $OUTDIR
#done
