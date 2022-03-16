###########Bench4############
BENCH4Q_BB=./dockq/only_backbone/bench4_dockqstats.csv
BENCH4Q_AA=./dockq/all_atoms/bench4_dockqstats.csv
BENCH4Q_RF=./dockq/only_backbone/bench4_dockqstats_RF.csv
BENCH4_KD=../../data/dockground/dockground.csv
#PCONSDOCK_BENCH4=./plddt/pconsdock-bench4.csv
PLDDT_BENCH4=./plddt/plddt_metrics_bench4_af_std_and_hhblits_msas_model_1_recycle_10_all_runs.csv
DSSP_BENCH4=../../data/dockground/dssp_df.csv
AFNEFF_BENCH4=../../data/dockground/AF_neffs.csv
TOPNEFF_BENCH4=../../data/dockground/top_neffs.csv
#############Marks positive#############
MARKSQ_RF=./dockq/only_backbone/dockqstats_marks_RF.csv
MARKSQ_AF_BB=./dockq/only_backbone/dockqstats_marks_af.csv
MARKSQ_AF_AA=./dockq/all_atoms/dockqstats_marks_af.csv
MARKSQ_GRAMM=./dockq/only_backbone/dockqstats_marks_gramm-aace18.csv #says only_backbone but is full atom
MARKSQ_TMFULL=./dockq/only_backbone/tmdockfull_dockqstats.csv
MARKSQ_TMINT=./dockq/only_backbone/tmdockint_dockqstats.csv
MARKSQ_MDOCKPP=./dockq/all_atoms/dockq_MdockPP.csv
PLDDT_MARKS_AF=./plddt/plddt_metrics_marks_af_std_and_hhblits_msas_model_1_recycle_10_all_runs.csv
#PCONSDOCK_MARKS=./plddt/pconsdock-marks.csv
PLDDT_MARKS_FUSED=./plddt/plddt_metrics_marks_af_std_and_fused_msas_model_1_recycle_10_all_runs.csv #Change this
DSSP_MARKS=../../data/marks/dssp_df.csv
IFSTATS_MARKS=../../data/marks/if_stats_marks.csv
ALN_SCORES_MARKS=../../data/marks/final_ids_lens.csv
OXSTATS_MARKS=../../data/marks/hhblits/oxstats.csv
AFNEFF_MARKS=../../data/marks/AF_neffs.csv
TOPNEFF_MARKS=../../data/marks/top_neffs.csv
AF_CHAIN_OVERLAP_MARKS=../../data/marks/overlap.csv
#############Marks negative#############
PLDDT_MARKS_AF_NEG=./plddt/plddt_metrics_marks_negative_fused_and_hhblits_msas.csv
#############Negatome#############
PLDDT_NEGATOME_AF_NEG=./plddt/plddt_metrics_negatome_fused_and_hhblits_msas.csv
#############New set#############
NEWSETQ=./dockq/only_backbone/dockqstats_newdimers_af.csv
PLDDT_NEWSET=./plddt/plddt_metrics_newset_af_std_and_hhblits_msas_model_1_recycle_10_all_runs.csv

OUTDIR=./plots/
python3 vis_analyze.py --bench4_dockq_aa $BENCH4Q_AA --bench4_dockq_RF $BENCH4Q_RF \
--plDDT_bench4 $PLDDT_BENCH4 --bench4_kingdom $BENCH4_KD --dssp_bench4 $DSSP_BENCH4 \
--afdefault_neff_bench4 $AFNEFF_BENCH4 --tophits_neff_bench4 $TOPNEFF_BENCH4 \
--marks_dockq_RF $MARKSQ_RF --marks_dockq_AF_bb $MARKSQ_AF_BB --marks_dockq_AF_aa $MARKSQ_AF_AA --marks_dockq_GRAMM $MARKSQ_GRAMM  \
--marks_dockq_TMfull $MARKSQ_TMFULL --marks_dockq_TMint $MARKSQ_TMINT \
--marks_dockq_mdockpp $MARKSQ_MDOCKPP \
--plDDT_marks_af $PLDDT_MARKS_AF --plDDT_marks_fused $PLDDT_MARKS_FUSED --dssp_marks $DSSP_MARKS \
--ifstats_marks $IFSTATS_MARKS --aln_scores_marks $ALN_SCORES_MARKS \
--oxstats_marks $OXSTATS_MARKS --afdefault_neff_marks $AFNEFF_MARKS --tophits_neff_marks $TOPNEFF_MARKS \
--af_chain_overlap_marks $AF_CHAIN_OVERLAP_MARKS \
--plDDT_marks_negative_af $PLDDT_MARKS_AF_NEG \
--plDDT_negatome_af $PLDDT_NEGATOME_AF_NEG \
--newset_dockq_AF $NEWSETQ --plDDT_newset $PLDDT_NEWSET \
--outdir $OUTDIR
