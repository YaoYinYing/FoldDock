# #Pairplots of hhblits vs AF msas
# convert DockQ_AF_vs_paired.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" DockQ_AF_vs_paired.svg
# convert DockQ_AF_vs_paired+AF2.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "C" DockQ_AF_vs_paired+AF2.svg
# montage DockQ_AF_vs_paired.svg DockQ_AF_vs_paired+AF2.svg -tile 2x1 -geometry +2+2 DockQ_merged_vs_AF.svg
# #Montage fig 1
# #Add labels
# convert DockQ_bench4.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" DockQ_bench4.svg
# montage DockQ_bench4.svg DockQ_merged_vs_AF.svg -tile 1x2 -geometry +2+2 Fig1.svg
# #Montage fig 2
# #Add labels
convert ROC_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" ROC_marks.svg
convert if_conctacts_vs_plddt.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" if_conctacts_vs_plddt.svg
convert pDockQ.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "C" pDockQ.svg
convert DockQ_marks_5runs.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "D" DockQ_marks_5runs.svg
# convert Fig2C_resized.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "C" Fig2C_labelled.svg
#
montage ROC_marks.svg if_conctacts_vs_plddt.svg \
pDockQ.svg DockQ_marks_5runs.svg -tile 2x2 -geometry +2+2 Fig2.svg

# #montage Fig2.svg Fig2C_labelled.svg -tile 1x2 -geometry +2+2 Fig2.svg
#
# #Montage msa comparison
# montage 5D1M_paired.svg 5D1M_af.svg 5D1M_both.svg -tile 3x1 -geometry +2+2 5D1M.svg
#
# #Montage fig 3
#Add labels
convert DockQ_per_SS_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" DockQ_per_SS_marks.svg
convert DockQ_per_num_if_contacts_total_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" DockQ_per_num_if_contacts_total_marks.svg
convert DockQ_per_tophit_Neff_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "C" DockQ_per_tophit_Neff_marks.svg
convert DockQ_per_org_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "D" DockQ_per_org_marks.svg
montage DockQ_per_SS_marks.svg DockQ_per_num_if_contacts_total_marks.svg  DockQ_per_tophit_Neff_marks.svg DockQ_per_org_marks.svg -tile 2x2 -geometry +2+2 Fig3.svg

#
# #Montage fig 4
# #Add labels
# convert RF_vs_AF_test.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" RF_vs_AF_test.svg
# convert GRAMM_vs_AF_test.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" GRAMM_vs_AF_test.svg
# montage RF_vs_AF_test.svg GRAMM_vs_AF_test.svg -tile 2x1 -geometry +2+2 Fig4.svg

#Figure 5
# convert 7EIV.svg -pointsize 40 -gravity NorthEast -annotate +0+0 "DockQ=0.76" 7EIV_l.svg
# convert 7EIV_l.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" 7EIV_l.svg
# convert 7MEZ.svg -pointsize 40 -gravity NorthEast -annotate +0+0 "DockQ=0.53" 7MEZ_l.svg
# convert 7MEZ_l.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" 7MEZ_l.svg
# convert 7EL1.svg -pointsize 40 -gravity NorthEast -annotate +0+0 "DockQ=0.01" 7EL1_l.svg
# convert 7EL1_l.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "C" 7EL1_l.svg
# convert 7LF7.svg -pointsize 40 -gravity NorthEast -annotate +0+0 "DockQ=0.02" 7LF7_l.svg
# convert 7LF7_l.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "D" 7LF7_l.svg
#
# montage 7EIV_l.svg 7MEZ_l.svg 7EL1_l.svg 7LF7_l.svg -tile 2x2 -geometry +2+2 Fig5.svg

#Box label
# convert DockQ_box_test.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" DockQ_box_test.svg
# #
# #Figure S1
#Kingdom stats
convert tophit_Neff_per_kd_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" tophit_Neff_per_kd_marks.svg
convert top_ranked_model_DockQ_af2_per_kd_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" top_ranked_model_DockQ_af2_per_kd_marks.svg
convert AFdefault_Neff_per_kd_marks.svg  -pointsize 100 -gravity NorthWest -annotate +0+0 "C" AFdefault_Neff_per_kd_marks.svg
montage  tophit_Neff_per_kd_marks.svg  top_ranked_model_DockQ_af2_per_kd_marks.svg AFdefault_Neff_per_kd_marks.svg -tile 3x1 -geometry +2+2 FigS1.svg

#Figure S2
#Add labels
convert DockQ_per_AFdefault_Neff_marks.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" DockQ_per_AFdefault_Neff_marks.svg
convert 'DockQ_per_biggest chain length_marks.svg' -pointsize 100 -gravity NorthWest -annotate +0+0 "B" 'DockQ_per_biggest chain length_marks.svg'
convert 'DockQ_per_smallest chain length_marks.svg' -pointsize 100 -gravity NorthWest -annotate +0+0 "C" 'DockQ_per_smallest chain length_marks.svg'
montage DockQ_per_AFdefault_Neff_marks.svg 'DockQ_per_biggest chain length_marks.svg' 'DockQ_per_smallest chain length_marks.svg'  -tile 3x1 -geometry +2+2 FigS2.svg

#Dev vs test
# convert dev_vs_test_max_chain_len.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" dev_vs_test_max_chain_len.svg
# convert dev_vs_test_min_chain_len.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" dev_vs_test_min_chain_len.svg
# convert dev_vs_test_num_if_contacts_total.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "C" dev_vs_test_num_if_contacts_total.svg
# convert dev_vs_test_AFdefault_Neff.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "D" dev_vs_test_AFdefault_Neff.svg
# convert dev_vs_test_tophit_Neff.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "E" dev_vs_test_tophit_Neff.svg
# montage dev_vs_test_max_chain_len.svg dev_vs_test_min_chain_len.svg dev_vs_test_num_if_contacts_total.svg \
#  dev_vs_test_AFdefault_Neff.svg dev_vs_test_tophit_Neff.svg -tile 3x2 -geometry +2+2 dev_vs_test.svg


#Marks pos vs negative
convert ROC_pos_neg.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" ROC_pos_neg.svg
montage if_plddt_av_pos_vs_neg_distr.svg  num_atoms_in_interface_pos_vs_neg_distr.svg \
pDockQ_pos_vs_neg_distr.svg -tile 2x2 -geometry +2+2 distr.svg
convert distr.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" distr.svg
montage ROC_pos_neg.svg distr.svg -tile 2x1 -geometry +2+2 Fig5.svg


#Figure S7 - Chain overlap
# convert overlap_distr.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "A" overlap_distr.svg
# convert dockq_vs_overlap.svg -pointsize 100 -gravity NorthWest -annotate +0+0 "B" dockq_vs_overlap.svg
# montage overlap_distr.svg dockq_vs_overlap.svg -tile 2x1 -geometry +2+2 FigS7.svg
