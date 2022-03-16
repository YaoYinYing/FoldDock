
#Evaluate the Docking
##########new_dimers#######
##AF2
DOCKQFILES=./dockqstats_newdimers_af2   #Path to file location
OUTFILE=./dockqstats_newdimers_af.csv
#python3 ./eval_docking.py --dockqfiles $DOCKQFILES --outfile $OUTFILE



##########Marks###########

##RF
DOCKQFILES=./dockqstats_marks_RF   #Path to file location
OUTFILE=./dockqstats_marks_RF.csv
#python3 ./eval_docking.py --dockqfiles $DOCKQFILES --outfile $OUTFILE

##AF2
DOCKQFILES=./dockqstats_marks_af2   #Path to file location
OUTFILE=./dockqstats_marks_af.csv
#python3 ./eval_docking.py --dockqfiles $DOCKQFILES --outfile $OUTFILE



#########Benchmark 4########
#Model 1
#Evaluate the docking of all different runs
DOCKQFILES=./dockqstats_bench4_af2   #Path to file location with dockq stats
OUTFILE=./bench4_dockqstats.csv #Where to write the merged csv with extracted dockq scores
#python3 ./eval_docking.py --dockqfiles $DOCKQFILES --outfile $OUTFILE

#Model 2
DOCKQFILES=./dockqstats_bench4_af2_af2stdmsa_model_2   #Path to file location with dockq stats
OUTFILE=./bench4_model2_dockqstats.csv #Where to write the merged csv with extracted dockq scores
python3 ./eval_docking.py --dockqfiles $DOCKQFILES --outfile $OUTFILE

#model_1, 10 recycles, both MSAs
# cd bench4_both_msas_model_1_rec10
# DOCKQFILES=./dockqstats_bench4   #Path to file location with dockq stats
# OUTFILE=./bench4_dockqstats_5runs.csv #Where to write the merged csv with extracted dockq scores
# python3 ../eval_docking.py --dockqfiles $DOCKQFILES --outfile $OUTFILE

#RoseTTAFold
DOCKQFILES=./dockqstats_bench4_RF   #Path to file location with dockq stats
OUTFILE=./bench4_dockqstats_RF.csv #Where to write the merged csv with extracted dockq scores
#python3 ./eval_docking.py --dockqfiles $DOCKQFILES --outfile $OUTFILE
