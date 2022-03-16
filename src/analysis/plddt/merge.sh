#bench4
PLDDTDIR=./plddt_metrics_bench4_af_std_and_hhblits_msas_model_1_recycle_10_
OUTNAME=./plddt_metrics_bench4_af_std_and_hhblits_msas_model_1_recycle_10_all_runs.csv
#python3 ./merge_plDDT_dfs.py --plDDTdir $PLDDTDIR --outname $OUTNAME

#Marks
#AF2+paired
PLDDTDIR=./plddt_metrics_marks_af_std_and_hhblits_msas_model_1_recycle_10_
OUTNAME=./plddt_metrics_marks_af_std_and_hhblits_msas_model_1_recycle_10_all_runs.csv
#python3 ./merge_plDDT_dfs.py --plDDTdir $PLDDTDIR --outname $OUTNAME

#Fused
PLDDTDIR=./plddt_metrics_marks_fused_and_hhblits_msas_model_1_recycle_10_
OUTNAME=./plddt_metrics_marks_af_std_and_fused_msas_model_1_recycle_10_all_runs.csv
python3 ./merge_plDDT_dfs.py --plDDTdir $PLDDTDIR --outname $OUTNAME

#New set
PLDDTDIR=./plddt_metrics_newset_af_std_and_hhblits
OUTNAME=./plddt_metrics_newset_af_std_and_hhblits_msas_model_1_recycle_10_all_runs.csv
#python3 ./merge_plDDT_dfs.py --plDDTdir $PLDDTDIR --outname $OUTNAME
